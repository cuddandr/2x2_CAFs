#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include/duneanaobj/StandardRecord/StandardRecord.h"

//Requires a file containing a list of input CAF files and returns an int code for success/error
//The list should be one file/path per line, and files/lines can be commented out using #
int test_plots(const std::string& file_list)
{
    std::vector<std::string> root_list;
    std::ifstream fin(file_list, std::ios::in);
    if(!fin.is_open())
    {
        std::cerr << "Failed to open " << file_list << std::endl;
        std::cerr << "Exiting" << std::endl;
        return 111;
    }
    else
    {
        std::cout << "Reading " << file_list << " for input ROOT files." << std::endl;
        std::string name;
        while(std::getline(fin, name))
        {
            if(name.front() == '#')
                continue;
            root_list.push_back(name);
        }
    }

    if(root_list.empty())
    {
        std::cerr << "No input ROOT files. Exiting." << std::endl;
        return 121;
    }

    TChain* caf_chain = new TChain("cafTree");
    for(const auto& file : root_list)
    {
        std::cout << "Adding " << file << " to TChain." << std::endl;
        caf_chain->Add(file.c_str());
    }

    //Beam direction -3.343 degrees in y
    const auto beam_dir = TVector3(0, -0.05836, 1.0);

    //Center of the 2x2 (if needed)
    const float tpc_x = 0.0;
    const float tpc_y = -268.0;
    const float tpc_z = (1333.5 + 1266.5) / 2.0;

    TH1F* h_num_ixn = new TH1F("h_num_ixn", "h_num_ixn", 10, 0, 10);

    TH1F* h_part_T = new TH1F("h_part_T", "h_part_T", 50, 0.0, 0.5);
    TH1F* h_part_T_cont = new TH1F("h_part_T_cont", "h_part_T_cont", 50, 0.0, 0.5);
    TH1F* h_part_p = new TH1F("h_part_p", "h_part_p", 50, 0.0, 1.0);
    TH1F* h_part_p_cont = new TH1F("h_part_p_cont", "h_part_p_cont", 50, 0.0, 1.0);
    TH1F* h_part_a = new TH1F("h_part_a", "h_part_a", 60, 0.0, 180);
    TH1F* h_part_a_cont = new TH1F("h_part_a_cont", "h_part_a_cont", 60, 0.0, 180);

    TH1F* h_vtx_x = new TH1F("h_vtx_x", "h_vtx_x", 50, -100 + tpc_x, 100 + tpc_x);
    TH1F* h_vtx_y = new TH1F("h_vtx_y", "h_vtx_y", 50, -100 + tpc_y, 100 + tpc_y);
    TH1F* h_vtx_z = new TH1F("h_vtx_z", "h_vtx_z", 50, -100 + tpc_z, 100 + tpc_z);

    auto sr = new caf::StandardRecord;
    caf_chain->SetBranchAddress("rec", &sr);

    const unsigned long nspills = caf_chain->GetEntries();
    const unsigned int incr = nspills / 10;
    std::cout << "Looping over " << nspills << " entries/spills..." << std::endl;

    const auto t_start{std::chrono::steady_clock::now()};
    for(unsigned long i = 0; i < nspills; ++i)
    {
        caf_chain->GetEntry(i);

        if(i % incr == 0)
            std::cout << "Spill #: " << i << std::endl;

        const auto num_ixn = sr->common.ixn.ndlp;
        h_num_ixn->Fill(num_ixn);

        //std::cout << "Num interactions: " << num_ixn << std::endl;
        for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
        {
            auto vtx = sr->common.ixn.dlp[ixn].vtx;

            if(std::isfinite(vtx.x) && std::isfinite(vtx.y) && std::isfinite(vtx.z))
            {
                h_vtx_x->Fill(vtx.x);
                h_vtx_y->Fill(vtx.y);
                h_vtx_z->Fill(vtx.z);
            }

            for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
            {
                auto part = sr->common.ixn.dlp[ixn].part.dlp[ipart];

                //if(!part.primary)
                //    continue;

                if(part.pdg == 13 or not part.primary)
                    continue;

                //if(std::abs(part.pdg) != 211)
                //    continue;

                float p = std::sqrt((part.p.x * part.p.x) + (part.p.y * part.p.y) + (part.p.z * part.p.z));

                auto dir = TVector3(part.end) - TVector3(part.start);
                auto angle = dir.Angle(beam_dir) * 180.0 / TMath::Pi();
                //auto angle = TVector3(part.p).Angle(beam_dir) * 180.0 / TMath::Pi();
                auto len = dir.Mag();

                h_part_T->Fill(part.E);
                h_part_p->Fill(p);
                h_part_a->Fill(angle);

                //if(part.contained)
                if(len > 2.0)
                {
                    h_part_T_cont->Fill(part.E);
                    h_part_p_cont->Fill(p);
                    h_part_a_cont->Fill(angle);
                }
            }
        }
    }
    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    TFile* caf_output = new TFile("caf_plots.root", "recreate");
    h_num_ixn->Write();
    h_vtx_x->Write();
    h_vtx_y->Write();
    h_vtx_z->Write();
    h_part_T->Write();
    h_part_T_cont->Write();
    h_part_p->Write();
    h_part_p_cont->Write();
    h_part_a->Write();
    h_part_a_cont->Write();
    caf_output->Close();

    std::cout << "Time elapsed: " << t_elapsed.count() << std::endl;
    std::cout << "Finished." << std::endl;
    return 0;
}
