#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include/duneanaobj/StandardRecord/StandardRecord.h"

int mesonless(const std::string& file_list)
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
    std::cout << "Processing " << root_list.size() << " files..." << std::endl;

    //Beam pointed 3.34 degrees in the -y direction (i.e. down)
    const auto beam_dir = TVector3(0.0, -0.05836, 1.0);

    //Offsets for 2x2 TPC center
    const float tpc_x = 0.0;
    const float tpc_y = -268.0;
    const float tpc_z = (1333.5 + 1266.5) / 2.0;

    TH1F* h_part_T = new TH1F("h_part_T", "h_part_T", 50, 0.0, 0.5);
    TH1F* h_part_T_cont = new TH1F("h_part_T_cont", "h_part_T_cont", 50, 0.0, 0.5);
    TH1F* h_part_p = new TH1F("h_part_p", "h_part_p", 50, 0.0, 1.0);
    TH1F* h_part_p_cont = new TH1F("h_part_p_cont", "h_part_p_cont", 50, 0.0, 1.0);
    TH1F* h_part_a = new TH1F("h_part_a", "h_part_a", 60, 0.0, 180);
    TH1F* h_part_a_cont = new TH1F("h_part_a_cont", "h_part_a_cont", 60, 0.0, 180);
    TH1F* h_part_l = new TH1F("h_part_l", "h_part_l", 50, 0.0, 100);
    TH1F* h_part_l_cont = new TH1F("h_part_l_cont", "h_part_l_cont", 50, 0.0, 100);

    TH1F* h_vtx_x = new TH1F("h_vtx_x", "h_vtx_x", 50, -100 + tpc_x, 100 + tpc_x);
    TH1F* h_vtx_y = new TH1F("h_vtx_y", "h_vtx_y", 50, -100 + tpc_y, 100 + tpc_y);
    TH1F* h_vtx_z = new TH1F("h_vtx_z", "h_vtx_z", 50, -100 + tpc_z, 100 + tpc_z);

    TH1F* h_num_muons = new TH1F("h_num_muons", "h_num_muons", 10, 0, 10);
    TH1F* h_num_pions = new TH1F("h_num_pions", "h_num_pions", 10, 0, 10);
    TH1F* h_num_prtns = new TH1F("h_num_prtns", "h_num_prtns", 10, 0, 10);
    TH1F* h_num_phtns = new TH1F("h_num_phtns", "h_num_phtns", 10, 0, 10);
    TH1F* h_num_elecs = new TH1F("h_num_elecs", "h_num_elecs", 10, 0, 10);
    TH1F* h_num_prims = new TH1F("h_num_prims", "h_num_prims", 10, 0, 10);

    auto sr = new caf::StandardRecord;
    caf_chain->SetBranchAddress("rec", &sr);

    unsigned int num_passed = 0;
    unsigned int total_ixn  = 0;
    unsigned int total_truth = 0;
    unsigned int total_rock  = 0;
    unsigned int total_argon = 0;

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
        for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
        {
            total_ixn++;
            const auto& vtx = sr->common.ixn.dlp[ixn].vtx;

            if(std::isfinite(vtx.x) && std::isfinite(vtx.y) && std::isfinite(vtx.z))
            {
                h_vtx_x->Fill(vtx.x);
                h_vtx_y->Fill(vtx.y);
                h_vtx_z->Fill(vtx.z);
            }

            const auto& vec_part = sr->common.ixn.dlp[ixn].part.dlp;
            unsigned int num_muons = 0;
            unsigned int num_prtns = 0;
            unsigned int num_pions = 0;
            unsigned int num_phtns = 0;
            unsigned int num_elecs = 0;
            unsigned int num_prims = 0;
            for(const auto& p : vec_part)
            {
                if(p.primary == false)
                    continue;
                num_prims++;

                if(p.pdg == 13)
                    num_muons++;

                if(p.pdg == 211)
                    num_pions++;

                if(p.pdg == 2212)
                    num_prtns++;

                if(p.pdg == 22)
                    num_phtns++;

                if(p.pdg == 11)
                    num_elecs++;
            }

            h_num_prims->Fill(num_prims);
            h_num_muons->Fill(num_muons);
            h_num_pions->Fill(num_pions);
            h_num_prtns->Fill(num_prtns);
            h_num_phtns->Fill(num_phtns);
            h_num_elecs->Fill(num_elecs);

            if((num_muons != 1) or (num_pions > 0))
                continue;

            if(num_prtns == 0)
                continue;

            num_passed++;
            std::cout << "Reco particles: ";
            std::cout << "[";
            for(const auto& p : vec_part)
            {
                if(p.primary == false)
                    continue;
                //std::cout << p.tgtA << "|";
                std::cout << p.pdg << ", ";
            }
            std::cout << "]\n";

            const auto& vec_truth_ixn = sr->common.ixn.dlp[ixn].truth;
            const auto& vec_overlap = sr->common.ixn.dlp[ixn].truthOverlap;
            for(unsigned long t = 0; t < vec_truth_ixn.size(); ++t)
            {
                auto idx = vec_truth_ixn.at(t);
                //std::cout << idx << std::endl;
                const auto& truth_ixn = sr->mc.nu[idx];
                std::cout << "True idx " << idx << "--> VtxID : " << truth_ixn.id << " w/ overlap " << vec_overlap.at(t) << std::endl;
                std::cout << "React : " << std::boolalpha << truth_ixn.iscc << "|" << truth_ixn.mode << std::endl;
                std::cout << "N prtn: " << truth_ixn.nproton << std::endl;
                std::cout << "N neut: " << truth_ixn.nneutron << std::endl;
                std::cout << "N pip : " << truth_ixn.npip << std::endl;
                std::cout << "N pim : " << truth_ixn.npim << std::endl;
                std::cout << "N pi0 : " << truth_ixn.npi0 << std::endl;

                std::cout << "Nprim : " << truth_ixn.nprim << std::endl;
                std::cout << "Nsec  : " << truth_ixn.nsec  << std::endl;

                if(truth_ixn.id >= 0 && truth_ixn.id < 1E9)
                    total_argon++;
                else
                    total_rock++;
                total_truth++;
            }

            for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
            {
                auto part = sr->common.ixn.dlp[ixn].part.dlp[ipart];

                //if((part.pdg != 2212 and part.pdg != 13) or not part.primary)
                if(part.pdg != 2212 or not part.primary)
                    continue;

                const auto& vec_truth_id = part.truth;
                const auto& vec_overlap  = part.truthOverlap;
                //for(const auto& id : vec_truth_id)
                std::cout << "Reco " << part.pdg << std::endl;
                for(unsigned int z = 0; z < vec_truth_id.size(); ++z)
                {
                    auto id = vec_truth_id.at(z);
                    caf::SRTrueParticle* q = nullptr;
                    if(id.type == 1)
                        q = &(sr->mc.nu[id.ixn].prim[id.part]);
                    if(id.type == 3)
                        q = &(sr->mc.nu[id.ixn].sec[id.part]);

                    if(q != nullptr)
                        std::cout << "-> " << q->pdg << " w/ overlap " << vec_overlap.at(z) << std::endl;
                }


                float p = std::sqrt((part.p.x * part.p.x) + (part.p.y * part.p.y) + (part.p.z * part.p.z));

                auto dir = TVector3(part.end) - TVector3(part.start);
                //auto angle = dir.Angle(beam_dir) * 180.0 / TMath::Pi();
                auto angle = TVector3(part.p).Angle(beam_dir) * 180.0 / TMath::Pi();
                auto len = dir.Mag();

                h_part_T->Fill(part.E);
                h_part_p->Fill(p);
                h_part_a->Fill(angle);
                h_part_l->Fill(len);

                if(part.contained && (len > 2.0))
                {
                    h_part_T_cont->Fill(part.E);
                    h_part_p_cont->Fill(p);
                    h_part_a_cont->Fill(angle);
                    h_part_l_cont->Fill(len);
                }
            }
        }
    }
    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    TFile* caf_output = new TFile("plots_mesonless.root", "recreate");
    h_num_prims->Write();
    h_num_muons->Write();
    h_num_pions->Write();
    h_num_prtns->Write();
    h_num_phtns->Write();
    h_num_elecs->Write();
    h_vtx_x->Write();
    h_vtx_y->Write();
    h_vtx_z->Write();
    h_part_T->Write();
    h_part_T_cont->Write();
    h_part_p->Write();
    h_part_p_cont->Write();
    h_part_a->Write();
    h_part_a_cont->Write();
    h_part_l->Write();
    h_part_l_cont->Write();
    caf_output->Close();

    float eff = (float)num_passed / (float)total_ixn;
    std::cout << "Selected " << num_passed << " / " << total_ixn << " ~ " << eff << std::endl;
    std::cout << "Truth total/rock/argon " << total_truth << "/" << total_rock << "/" << total_argon << std::endl;
    std::cout << "Time elapsed: " << t_elapsed.count() << std::endl;
    std::cout << "Finished." << std::endl;
    return 0;
}
