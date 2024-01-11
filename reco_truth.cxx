#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include/duneanaobj/StandardRecord/StandardRecord.h"

int reco_truth(const std::string& file_list)
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
            const auto& vec_truth_ixn = sr->common.ixn.dlp[ixn].truth;
            const auto& vec_overlap = sr->common.ixn.dlp[ixn].truthOverlap;

            std::cout << "----------------------------------------" << std::endl;
            std::cout << "Spill " << i << "; Reco ixn " << ixn << std::endl;
            bool is_cc0pi = false;
            for(unsigned long t = 0; t < vec_truth_ixn.size(); ++t)
            {
                auto idx = vec_truth_ixn.at(t);
                const auto& truth_ixn = sr->mc.nu[idx];

                const unsigned int npi = truth_ixn.npip + truth_ixn.npim + truth_ixn.npi0;
                if(truth_ixn.nproton >= 0 && npi == 0 && truth_ixn.iscc)
                {
                    is_cc0pi = true;
                    num_passed += 1;

                    std::cout << "True idx " << idx << std::endl;
                    std::cout << "VtxID : " << truth_ixn.id << " w/ overlap " << vec_overlap.at(t) << std::endl;
                    std::cout << "isRock: " << std::boolalpha << (truth_ixn.id < 1E9 ? false : true) << std::endl;
                    std::cout << "isCC  : " << std::boolalpha << truth_ixn.iscc << std::endl;
                    std::cout << "Reac #: " << truth_ixn.mode << std::endl;
                    std::cout << "N prtn: " << truth_ixn.nproton << std::endl;
                    std::cout << "N neut: " << truth_ixn.nneutron << std::endl;
                    std::cout << "N pip : " << truth_ixn.npip << std::endl;
                    std::cout << "N pim : " << truth_ixn.npim << std::endl;
                    std::cout << "N pi0 : " << truth_ixn.npi0 << std::endl;
                    std::cout << "Nprim : " << truth_ixn.nprim << std::endl;
                    std::cout << "Nsec  : " << truth_ixn.nsec  << std::endl;

                    std::cout << "Truth Primaries: " << std::endl;
                    std::cout << "\tG4ID | PDG | P_mag" << std::endl;
                    for(unsigned int p = 0; p < truth_ixn.nprim; ++p)
                    {
                        auto q = truth_ixn.prim[p];
                        std::cout << "\t" << q.G4ID << " | " << q.pdg << " | " << q.p.Mag() << std::endl;
                    }

                    std::cout << "Truth Secondaries: " << std::endl;
                    std::cout << "\tG4ID | PDG | P_mag" << std::endl;
                    for(unsigned int s = 0; s < truth_ixn.nsec; ++s)
                    {
                        auto q = truth_ixn.sec[s];
                        std::cout << "\t" << q.G4ID << " | " << q.pdg << " | " << q.p.Mag() << std::endl;
                    }

                }
            }

            if(!is_cc0pi)
                continue;

            std::cout << "Found CC0pi truth interaction" << std::endl;
            std::cout << "Printing selected/reco information" << std::endl;
            const auto& vtx = sr->common.ixn.dlp[ixn].vtx;
            const auto& vec_part = sr->common.ixn.dlp[ixn].part.dlp;

            /*
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
            */

            //if((num_muons != 1) or (num_pions > 0))
            //   continue;

            //if(num_prtns == 0)
            //    continue;

            std::cout << "Reco primaries (PDGs): ";
            std::cout << "[";
            for(const auto& p : vec_part)
            {
                if(p.primary == false)
                    continue;
                //std::cout << p.tgtA << "|";
                std::cout << p.pdg << ", ";
            }
            std::cout << "]\n";

            std::cout << "Reco secondaries (PDGs): ";
            std::cout << "[";
            for(const auto& p : vec_part)
            {
                if(p.primary == true)
                    continue;
                //std::cout << p.tgtA << "|";
                std::cout << p.pdg << ", ";
            }
            std::cout << "]\n";

            std::cout << "Printing reco/truth particle overlap" << std::endl;
            for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
            {
                auto part = sr->common.ixn.dlp[ixn].part.dlp[ipart];

                const auto& vec_truth_id = part.truth;
                const auto& vec_overlap  = part.truthOverlap;

                std::cout << "Reco " << ipart << " | " << (part.primary ? "Prim" : "Sec") << " | " << part.pdg << std::endl;
                /* std::cout << "Truth particle overlaps\n"; */
                std::cout << "\tG4ID | Type | PDG | Overlap" << std::endl;
                for(unsigned int z = 0; z < vec_truth_id.size(); ++z)
                {
                    auto id = vec_truth_id.at(z);
                    caf::SRTrueParticle* q = nullptr;
                    if(id.type == 1)
                        q = &(sr->mc.nu[id.ixn].prim[id.part]);
                    if(id.type == 3)
                        q = &(sr->mc.nu[id.ixn].sec[id.part]);

                    if(q != nullptr)
                        std::cout << "\t" << q->G4ID << " | " << id.type << " | " << q->pdg << " | " << vec_overlap.at(z) << std::endl;
                }
            }
        }
    }
    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    float eff = (float)num_passed / (float)total_ixn;
    std::cout << "Selected " << num_passed << " / " << total_ixn << " ~ " << eff << std::endl;
    //std::cout << "Truth total/rock/argon " << total_truth << "/" << total_rock << "/" << total_argon << std::endl;
    std::cout << "Time elapsed: " << t_elapsed.count() << std::endl;
    std::cout << "Finished." << std::endl;
    return 0;
}
