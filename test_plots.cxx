#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_01_00/include/duneanaobj/StandardRecord/StandardRecord.h"

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

    const float tpc_x = 0.0;
    const float tpc_y = -268.0;
    const float tpc_z = (1333.5 + 1266.5) / 2.0;

    TH1F* h_part_T = new TH1F("h_part_T", "h_part_T", 50, 0.0, 0.5);
    TH1F* h_part_T_cont = new TH1F("h_part_T_cont", "h_part_T_cont", 50, 0.0, 0.5);
    TH1F* h_part_p = new TH1F("h_part_p", "h_part_p", 50, 0.0, 1.0);
    TH1F* h_part_p_cont = new TH1F("h_part_p_cont", "h_part_p_cont", 50, 0.0, 1.0);

    TH1F* h_vtx_x = new TH1F("h_vtx_x", "h_vtx_x", 50, -100 + tpc_x, 100 + tpc_x);
    TH1F* h_vtx_y = new TH1F("h_vtx_y", "h_vtx_y", 50, -100 + tpc_y, 100 + tpc_y);
    TH1F* h_vtx_z = new TH1F("h_vtx_z", "h_vtx_z", 50, -100 + tpc_z, 100 + tpc_z);

    auto sr = new caf::StandardRecord;
    caf_chain->SetBranchAddress("rec", &sr);

    const unsigned long nspills = caf_chain->GetEntries();
    const unsigned int incr = nspills / 10;
    std::cout << "Looping over " << nspills << " entries/spills..." << std::endl;
    for(unsigned long i = 0; i < nspills; ++i)
    {
        caf_chain->GetEntry(i);

        if(i % incr == 0)
            std::cout << "Spill #: " << i << std::endl;

        //std::cout << "Spill #: " << i << std::endl;
        //std::cout << "Num interactions: " << sr->common.ixn.ndlp << std::endl;
        for(unsigned long ixn = 0; ixn < sr->common.ixn.ndlp; ++ixn)
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

                if(part.pdg != 2212)
                    continue;

                //if(std::abs(part.pdg) != 211)
                //    continue;

                float p = std::sqrt((part.p.x * part.p.x) + (part.p.y * part.p.y) + (part.p.z * part.p.z));

                h_part_T->Fill(part.E);
                h_part_p->Fill(p);

                if(part.contained)
                {
                    h_part_T_cont->Fill(part.E);
                    h_part_p_cont->Fill(p);
                }
            }
        }
    }

    TFile* caf_output = new TFile("caf_plots.root", "recreate");
    h_vtx_x->Write();
    h_vtx_y->Write();
    h_vtx_z->Write();
    h_part_T->Write();
    h_part_T_cont->Write();
    h_part_p->Write();
    h_part_p_cont->Write();
    caf_output->Close();

    std::cout << "Finished." << std::endl;
    return 0;
}
