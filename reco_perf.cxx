#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include/duneanaobj/StandardRecord/StandardRecord.h"

//Requires a file containing a list of input CAF files and returns an int code for success/error
//The list should be one file/path per line, and files/lines can be commented out using #
int reco_perf(const std::string& file_list)
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
    std::cout << "Finished adding files..." << std::endl;

    //Beam direction -3.343 degrees in y
    const auto beam_dir = TVector3(0, -0.05836, 1.0);

    //Center of the 2x2 LAr
    //For MR4.5 CAFs the 2x2 coordinates are actually centered somewhere else
    const float tpc_x = 0.0;
    const float tpc_y = 0.0; // -268.0;
    const float tpc_z = 0.0; //(1333.5 + 1266.5) / 2.0;

    TH1F* h_num_ixn = new TH1F("h_num_ixn", "h_num_ixn", 10, 0, 10);

    TH1F* h_part_T = new TH1F("h_part_T", "h_part_T", 50, 0.0, 0.5);
    TH1F* h_part_p = new TH1F("h_part_p", "h_part_p", 50, 0.0, 1.0);
    TH1F* h_part_a = new TH1F("h_part_a", "h_part_a", 60, 0.0, 180);
    TH1F* h_true_T = new TH1F("h_true_T", "h_true_T", 50, 0.0, 0.5);
    TH1F* h_true_p = new TH1F("h_true_p", "h_true_p", 50, 0.0, 1.0);
    TH1F* h_true_a = new TH1F("h_true_a", "h_true_a", 60, 0.0, 180);

    TH1F* h_vtx_x = new TH1F("h_vtx_x", "h_vtx_x", 50, -100 + tpc_x, 100 + tpc_x);
    TH1F* h_vtx_y = new TH1F("h_vtx_y", "h_vtx_y", 50, -100 + tpc_y, 100 + tpc_y);
    TH1F* h_vtx_z = new TH1F("h_vtx_z", "h_vtx_z", 50, -100 + tpc_z, 100 + tpc_z);

    //Attach SR object to the tree, currently not using SRProxy (will swtich eventually)
    auto sr = new caf::StandardRecord;
    caf_chain->SetBranchAddress("rec", &sr);

    const unsigned long nspills = caf_chain->GetEntries();
    const unsigned int incr = nspills / 10;
    std::cout << "Looping over " << nspills << " entries/spills..." << std::endl;

    //First loop over each spill, then each reco interaction, then each reco particle
    const auto t_start{std::chrono::steady_clock::now()};
    for(unsigned long i = 0; i < nspills; ++i)
    {
        caf_chain->GetEntry(i);

        if(i % incr == 0)
            std::cout << "Spill #: " << i << std::endl;

        const auto num_ixn = sr->common.ixn.ndlp;
        h_num_ixn->Fill(num_ixn);

        for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
        {
            auto vtx = sr->common.ixn.dlp[ixn].vtx;

            //Get the truth interaction(s) corresponding to this reco interaction
            const auto& vec_truth_ixn = sr->common.ixn.dlp[ixn].truth;
            const auto& vec_overlap_ixn = sr->common.ixn.dlp[ixn].truthOverlap;

            if(vec_overlap_ixn.empty())
                continue;

            //Find the truth interaction with the largest overlap
            auto result = std::max_element(vec_overlap_ixn.begin(), vec_overlap_ixn.end());
            auto max_overlap = std::distance(vec_overlap_ixn.begin(), result);
            auto truth_idx = vec_truth_ixn.at(max_overlap);
            auto truth_ixn = sr->mc.nu[truth_idx];

            //Put cuts on true interaction quantities here
            //For example, reject interactions not on argon
            if(truth_ixn.targetPDG != 1000180400)
                continue;

            if(std::isfinite(vtx.x) && std::isfinite(vtx.y) && std::isfinite(vtx.z))
            {
                h_vtx_x->Fill(vtx.x);
                h_vtx_y->Fill(vtx.y);
                h_vtx_z->Fill(vtx.z);
            }

            for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
            {
                //Store current reco particle for easier access
                auto part = sr->common.ixn.dlp[ixn].part.dlp[ipart];

                //Put reco particle cuts here
                if(part.pdg == 13 or not part.primary)
                    continue;

                //Get truth match(es) for this reco particle
                caf::SRTrueParticle* truth_match = nullptr;
                const auto& vec_truth_id = part.truth;
                const auto& vec_overlap  = part.truthOverlap;

                //If the truth overlap vector is empty, then assume no truth match and skip
                if(vec_overlap.empty())
                {
                    //std::cout << "No truth match... skipping reco particle..." << std::endl;
                    continue;
                }

                //Find the truth particle with the largest overlap
                auto result = std::max_element(vec_overlap.begin(), vec_overlap.end());
                auto max_overlap = std::distance(vec_overlap.begin(), result);
                auto truth_id = vec_truth_id.at(max_overlap);

                //Get pointer to the corresponding truth particle
                if(truth_id.type == 1)
                    truth_match = &(sr->mc.nu[truth_id.ixn].prim[truth_id.part]);
                else if(truth_id.type == 3)
                    truth_match = &(sr->mc.nu[truth_id.ixn].sec[truth_id.part]);
                else
                {
                    std::cout << "Invalid truth id type!" << std::endl;
                    continue;
                }

                //Finally get or calculate various reco/truth quantities
                float p = std::sqrt((part.p.x * part.p.x) + (part.p.y * part.p.y) + (part.p.z * part.p.z));
                auto dir = TVector3(part.end) - TVector3(part.start);
                auto angle = dir.Angle(beam_dir) * 180.0 / TMath::Pi(); //calculate angle w.r.t neutrino beam direction
                //auto angle = TVector3(part.p).Angle(beam_dir) * 180.0 / TMath::Pi(); //different method for the angle

                auto true_dir = TVector3(truth_match->end_pos) - TVector3(truth_match->start_pos);
                auto true_angle = true_dir.Angle(beam_dir) * 180.0 / TMath::Pi();

                //h_part_T->Fill(part.E);
                h_part_p->Fill(part.p.Mag());
                h_part_a->Fill(angle);

                //h_true_T->Fill(truth_match->E);
                h_true_p->Fill(truth_match->p.Mag());
                h_true_a->Fill(true_angle);
            }
        }
    }
    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    //Write plots to file, overwriting the ROOT file
    TFile* caf_output = new TFile("caf_plots.root", "recreate");
    h_num_ixn->Write();
    h_vtx_x->Write();
    h_vtx_y->Write();
    h_vtx_z->Write();
    h_part_T->Write();
    h_part_p->Write();
    h_part_a->Write();
    h_true_T->Write();
    h_true_p->Write();
    h_true_a->Write();
    caf_output->Close();

    std::cout << "Time elapsed: " << t_elapsed.count() << std::endl;
    std::cout << "Finished." << std::endl;
    return 0;
}
