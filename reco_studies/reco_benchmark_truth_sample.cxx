/* #include "/cvmfs/dune.opensciencegrid.org/products/dune/srproxy/v00.43/include/SRProxy/BasicTypesProxy.h" */
#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include/duneanaobj/StandardRecord/StandardRecord.h"
#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include/duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include <TLatex.h>
#include <TCanvas.h>
//Requires a file containing a list of input CAF files and returns an int code for success/error
//The list should be one file/path per line, and files/lines can be commented out using #
//Boolean flag to switch behavior for reading structured/flat CAFs (defaults to structured CAFs)
int reco_benchmark_truth_sample(const std::string& file_list, bool is_flat = false)
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

    // DEFINE: Vectors to hold information to keep in output TTree file 
    std::vector< double >  true_energy;
    std::vector< double >  true_p_x; 
    std::vector< double >  true_p_y; 
    std::vector< double >  true_p_z;
    std::vector< double >  true_p_mag;
    std::vector< double >  true_length;
    std::vector< double >  true_angle;
    std::vector< double >  true_angle_rot;
    std::vector< double >  true_angle_incl;
    std::vector< double >  true_angle_x;
    std::vector< double >  true_angle_y;
    std::vector< double >  true_angle_z;
    std::vector< double >  true_track_start_x;
    std::vector< double >  true_track_start_y;
    std::vector< double >  true_track_start_z;
    std::vector< double >  true_track_end_x;
    std::vector< double >  true_track_end_y;
    std::vector< double >  true_track_end_z;
    std::vector< int >     true_pdg;
    std::vector< int >     true_ixn_charged_track_mult;

    std::vector< int >     spill_index;
    std::vector< int >     file_index;
    std::vector< int >     event;
    std::vector< int >     run;
    std::vector< int >     subrun;
    std::vector< std::string > caf_file_name;
    

    // DEFINE: TTree and TBranches to go in output ROOT file
    TTree *fRecoBenchmarkTree=new TTree("RecoBenchmarkTree", "Reco Benchmark Variables");
    fRecoBenchmarkTree->Branch("true_energy", &true_energy);
    fRecoBenchmarkTree->Branch("true_p_x", &true_p_x);
    fRecoBenchmarkTree->Branch("true_p_y", &true_p_y);
    fRecoBenchmarkTree->Branch("true_p_z", &true_p_z);
    fRecoBenchmarkTree->Branch("true_p_mag", &true_p_mag);
    fRecoBenchmarkTree->Branch("true_length", &true_length);
    fRecoBenchmarkTree->Branch("true_angle", &true_angle);
    fRecoBenchmarkTree->Branch("true_angle_rot", &true_angle_rot);
    fRecoBenchmarkTree->Branch("true_angle_incl", &true_angle_incl);
    fRecoBenchmarkTree->Branch("true_angle_x", &true_angle_x);
    fRecoBenchmarkTree->Branch("true_angle_y", &true_angle_y);
    fRecoBenchmarkTree->Branch("true_angle_z", &true_angle_z);
    fRecoBenchmarkTree->Branch("true_track_start_x", &true_track_start_x);
    fRecoBenchmarkTree->Branch("true_track_start_y", &true_track_start_y);
    fRecoBenchmarkTree->Branch("true_track_start_z", &true_track_start_z);
    fRecoBenchmarkTree->Branch("true_track_end_x", &true_track_end_x);
    fRecoBenchmarkTree->Branch("true_track_end_y", &true_track_end_y);
    fRecoBenchmarkTree->Branch("true_track_end_z", &true_track_end_z);
    fRecoBenchmarkTree->Branch("true_pdg", &true_pdg);
    fRecoBenchmarkTree->Branch("true_ixn_charged_track_mult", &true_ixn_charged_track_mult);

    fRecoBenchmarkTree->Branch("event", &event);
    fRecoBenchmarkTree->Branch("run", &run);
    fRecoBenchmarkTree->Branch("subrun", &subrun);
    fRecoBenchmarkTree->Branch("caf_file_name", &caf_file_name);



    //Beam direction -3.343 degrees in y
    const auto beam_dir = TVector3(0, -0.05836, 1.0);

    //z-direction (roughly beam dir)
    const auto z_plus_dir = TVector3(0, 0, 1.0);
    const auto y_plus_dir = TVector3(0, 1.0, 0.0);
    const auto x_plus_dir = TVector3(1.0, 0, 0.0);

    //negative y-direction 
    const auto y_minus_dir = TVector3(0, -1.0, 0.0);

    //Center of the 2x2 LAr
    //For MR4.5 CAFs the 2x2 coordinates are actually centered somewhere else
    const float tpc_x = 0.0;
    const float tpc_y = 0.0; // -268.0;
    const float tpc_z = 0.0; //(1333.5 + 1266.5) / 2.0;

    
    //Attach SR object to the tree, not using SRProxy (SR Proxy used below)
    //auto sr = new caf::StandardRecord;
    //caf_chain->SetBranchAddress("rec", &sr);

    //const unsigned long nspills = caf_chain->GetEntries();
    //const unsigned int incr = nspills / 10;
    //std::cout << "Looping over " << nspills << " entries/spills..." << std::endl;

    //Most outer loop over files for SRProxy
    const auto t_start{std::chrono::steady_clock::now()};
    auto file_num = 0;
    for(const auto& f : root_list)
    {
    std::cout << "Processing " << f << std::endl;
    file_num++;

    //Open file and attach SRProxy object
    //Different tree name for structured and flat CAFs
    TFile* caf_file = TFile::Open(f.c_str(), "READ");
    TTree* caf_tree = (TTree*)caf_file->Get("cafTree");
    std::string tree_name = is_flat ? "rec" : "";
    auto sr = new caf::SRProxy(caf_tree, tree_name);

    //First loop over each spill, then each reco interaction, then each reco particle
    const unsigned long nspills = caf_tree->GetEntries();
    const unsigned int incr = nspills / 10;
    std::cout << "Looping over " << nspills << " entries/spills..." << std::endl;
    for(unsigned long i = 0; i < nspills; ++i)
    {
        caf_tree->GetEntry(i);

        auto spill_num = i;

        if(i % incr == 0)
            std::cout << "Spill #: " << i << std::endl;

        const auto num_ixn = sr->mc.nu.size();

        for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
        {

            const auto& truth_ixn = sr->mc.nu[ixn];

            //Put cuts on true interaction quantities here
            //For example, reject interactions not on argon
            if(truth_ixn.targetPDG != 1000180400)
                continue;

            //if(std::isfinite(vtx.x) && std::isfinite(vtx.y) && std::isfinite(vtx.z))
            //{
            //    //enter action here
            //}

            auto true_ixn_trks = 0;
            // Loop over particles in matched true interaction
            for(unsigned long ipart = 0; ipart < sr->mc.nu[ixn].prim.size(); ++ipart)
            {
                //Store current true particle for easier access
                const auto& true_part = sr->mc.nu[ixn].prim[ipart];
                if((true_part.pdg != 2212 and true_part.pdg != 321 and true_part.pdg != 211 and true_part.pdg != 13))
                    continue;

                //Put true particle cuts here
                if((abs(abs(true_part.start_pos.z)-64.538) < 1.0) and (abs(abs(true_part.end_pos.z)-64.538)) < 1.0){
                    //std::cout<<"Particle Start:"<<part.start.z<<" End:"<<part.end.z<<std::endl;
                    continue;
                }

                //Put true particle cuts here
                if((abs(true_part.start_pos.z+64.538) < 1.0) or (abs(true_part.end_pos.z+64.538) < 1.0)){
                    //std::cout<<"Particle Start:"<<part.start.z<<" End:"<<part.end.z<<std::endl;
                    continue;
                }

                ++true_ixn_trks;

                //Finally get or calculate various truth quantities

                auto true_pvec = TVector3(true_part.p.px, true_part.p.py, true_part.p.pz);
                auto true_dir = TVector3(true_part.end_pos.x, true_part.end_pos.y, true_part.end_pos.z)
                                - TVector3(true_part.start_pos.x, true_part.start_pos.y, true_part.start_pos.z);
                std::cout<<"True End X: "<<true_part.end_pos.x<<std::endl; 
                auto true_cos_angle = TMath::Cos(true_dir.Angle(beam_dir));
                true_dir.RotateY(-TMath::Pi()/2);
                auto true_cos_rot_anode_angle = TMath::Cos(true_dir.Theta()); //calculate cosine of track rotational angle (projection on anode)
                auto true_cos_incl_anode_angle = TMath::Cos(true_dir.Phi()); //calculate cosine of track inclination angle (off of anode)
                true_dir.RotateY(TMath::Pi()/2);
                auto true_length_val = true_dir.Mag();

                true_dir.RotateY(-TMath::Pi()/2);

                // POPULATE: Record information in vectors (defined above) for tracks that
	            //           have passed all cuts

                true_energy.push_back(true_part.p.E);
                true_p_x.push_back(true_part.p.px); 
                true_p_y.push_back(true_part.p.py); 
                true_p_z.push_back(true_part.p.pz);
                true_p_mag.push_back(true_pvec.Mag());
                true_length.push_back(true_length_val);
                true_dir.RotateY(TMath::Pi()/2);
                true_angle.push_back(true_dir.Angle(beam_dir));
                true_angle_x.push_back(true_dir.Angle(x_plus_dir));
                true_angle_y.push_back(true_dir.Angle(y_plus_dir));
                true_angle_z.push_back(true_dir.Angle(z_plus_dir));
                true_dir.RotateY(-TMath::Pi()/2);
                true_angle_rot.push_back(true_dir.Theta());
                true_angle_incl.push_back(true_dir.Phi());
                true_track_start_x.push_back(true_part.start_pos.x);
                true_track_start_y.push_back(true_part.start_pos.y);
                true_track_start_z.push_back(true_part.start_pos.z);
                true_track_end_x.push_back(true_part.end_pos.x);
                true_track_end_y.push_back(true_part.end_pos.y);
                true_track_end_z.push_back(true_part.end_pos.z);
                true_pdg.push_back(true_part.pdg);
                
                event.push_back(sr->meta.nd_lar.event);
                run.push_back(sr->meta.nd_lar.run);
                subrun.push_back(sr->meta.nd_lar.subrun);
                caf_file_name.push_back(f.c_str());


            }
            // Loop over particles in the interaction again to load charged track multiplicity
            for(unsigned long ipart = 0; ipart < true_ixn_trks; ++ipart)
            {
                true_ixn_charged_track_mult.push_back(true_ixn_trks);
            }
        }
    } //end of reco interactions

    } //end of file loop
    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

        // Output TTree file name
    std::string file_name = "track_reco_benchmark_truth_sample";

    // DEFINE: Output TFile
    TFile *f=new TFile(Form("%s.root", file_name.c_str()),"RECREATE");

    // POPULATE: Fill TTree and write to output ROOT file
    fRecoBenchmarkTree->Fill();
    fRecoBenchmarkTree->Write();
    
    std::cout << "Filled and wrote TTree." << std::endl;

    // CLOSE: Output ROOT file
    f->Close();

    std::cout << "Time elapsed: " << t_elapsed.count() << std::endl;
    std::cout << "Finished." << std::endl;
    return 0;
}
