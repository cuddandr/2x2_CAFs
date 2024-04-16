/* #include "/cvmfs/dune.opensciencegrid.org/products/dune/srproxy/v00.43/include/SRProxy/BasicTypesProxy.h" */
#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include/duneanaobj/StandardRecord/StandardRecord.h"
#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include/duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include <TLatex.h>
#include <TCanvas.h>
//Requires a file containing a list of input CAF files and returns an int code for success/error
//The list should be one file/path per line, and files/lines can be commented out using #
//Boolean flag to switch behavior for reading structured/flat CAFs (defaults to structured CAFs)
int reco_benchmark(const std::string& file_list, bool is_flat = false)
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
    std::vector< double >  reco_energy;
    std::vector< double >  reco_p_x; 
    std::vector< double >  reco_p_y; 
    std::vector< double >  reco_p_z;
    std::vector< double >  reco_p_mag;
    std::vector< double >  reco_length;
    std::vector< double >  reco_angle;
    std::vector< double >  reco_angle_rot;
    std::vector< double >  reco_angle_incl;
    std::vector< double >  reco_angle_x;
    std::vector< double >  reco_angle_y;
    std::vector< double >  reco_angle_z;
    std::vector< double >  reco_track_start_x;
    std::vector< double >  reco_track_start_y;
    std::vector< double >  reco_track_start_z;
    std::vector< double >  reco_track_end_x;
    std::vector< double >  reco_track_end_y;
    std::vector< double >  reco_track_end_z;
    std::vector< int >     reco_pdg;
    std::vector< int >     reco_ixn_charged_track_mult;

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

    std::vector< double >  overlap;
    std::vector< double >  true_ixn_index;
    std::vector< double >  reco_ixn_index;
    std::vector< int >     spill_index;
    std::vector< int >     file_index;
    

    // DEFINE: TTree and TBranches to go in output ROOT file
    TTree *fRecoBenchmarkTree=new TTree("RecoBenchmarkTree", "Reco Benchmark Variables");
    fRecoBenchmarkTree->Branch("reco_energy", &reco_energy);
    fRecoBenchmarkTree->Branch("reco_p_x", &reco_p_x);
    fRecoBenchmarkTree->Branch("reco_p_y", &reco_p_y);
    fRecoBenchmarkTree->Branch("reco_p_z", &reco_p_z);
    fRecoBenchmarkTree->Branch("reco_p_mag", &reco_p_mag);
    fRecoBenchmarkTree->Branch("reco_length", &reco_length);
    fRecoBenchmarkTree->Branch("reco_angle", &reco_angle);
    fRecoBenchmarkTree->Branch("reco_angle_rot", &reco_angle_rot);
    fRecoBenchmarkTree->Branch("reco_angle_incl", &reco_angle_incl);
    fRecoBenchmarkTree->Branch("reco_angle_x", &reco_angle_x);
    fRecoBenchmarkTree->Branch("reco_angle_y", &reco_angle_y);
    fRecoBenchmarkTree->Branch("reco_angle_z", &reco_angle_z);
    fRecoBenchmarkTree->Branch("reco_track_start_x", &reco_track_start_x);
    fRecoBenchmarkTree->Branch("reco_track_start_y", &reco_track_start_y);
    fRecoBenchmarkTree->Branch("reco_track_start_z", &reco_track_start_z);
    fRecoBenchmarkTree->Branch("reco_track_end_x", &reco_track_end_x);
    fRecoBenchmarkTree->Branch("reco_track_end_y", &reco_track_end_y);
    fRecoBenchmarkTree->Branch("reco_track_end_z", &reco_track_end_z);
    fRecoBenchmarkTree->Branch("reco_pdg", &reco_pdg);
    fRecoBenchmarkTree->Branch("reco_ixn_charged_track_mult", &reco_ixn_charged_track_mult);
    fRecoBenchmarkTree->Branch("reco_ixn_index", &reco_ixn_index);


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
    fRecoBenchmarkTree->Branch("true_ixn_index", &true_ixn_index);

    fRecoBenchmarkTree->Branch("spill_index", &spill_index);
    fRecoBenchmarkTree->Branch("file_index", &file_index);

    fRecoBenchmarkTree->Branch("overlap", &overlap);


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

    TH1F* h_num_ixn = new TH1F("h_num_ixn", "h_num_ixn", 10, 0, 10);

    TH1F* h_part_T = new TH1F("h_part_T", "h_part_T", 50, 0.0, 0.5);
    TH1F* h_part_p = new TH1F("h_part_p", "h_part_p", 50, 0.0, 1.0);
    TH1F* h_part_a = new TH1F("h_part_a", "h_part_a", 50, -1.0,1.0);
    TH1F* h_part_a_rot = new TH1F("h_part_a_rot", "h_part_a_rot", 50, -1.0,1.0);
    TH1F* h_part_a_incl = new TH1F("h_part_a_incl", "h_part_a_incl", 50, -1.0,1.0);
    TH1F* h_part_l = new TH1F("h_part_l", "h_part_l", 30, 0.0, 60);
    TH1F* h_part_pdg = new TH1F("h_part_pdg", "h_part_pdg", 20000, -10000, 10000);

    TH1F* h_true_T = new TH1F("h_true_T", "h_true_T", 50, 0.0, 0.5);
    TH1F* h_true_p = new TH1F("h_true_p", "h_true_p", 50, 0.0, 1.0);
    TH1F* h_true_a = new TH1F("h_true_a", "h_true_a", 50, -1.0,1.0);
    TH1F* h_true_a_rot = new TH1F("h_true_a_rot", "h_true_a_rot", 50, -1.0,1.0);
    TH1F* h_true_a_incl = new TH1F("h_true_a_incl", "h_true_a_incl", 50, -1.0,1.0);
    TH1F* h_true_l = new TH1F("h_true_l", "h_true_l", 30, 0.0, 60);
    TH1F* h_true_pdg = new TH1F("h_true_pdg", "h_true_pdg", 20000, -10000, 10000);

    TH1F* h_diff_T = new TH1F("h_diff_T", "h_diff_T", 150, -0.5, 1.0);
    TH1F* h_diff_p = new TH1F("h_diff_p", "h_diff_p", 150, -1.0, 2.0);
    TH1F* h_diff_a = new TH1F("h_diff_a", "h_diff_a", 100, -2.0,2.0);
    TH1F* h_diff_a_rot = new TH1F("h_diff_a_rot", "h_diff_a_rot", 100, -2.0,2.0);
    TH1F* h_diff_a_incl = new TH1F("h_diff_a_incl", "h_diff_a_incl", 100, -2.0,2.0);
    TH1F* h_diff_l = new TH1F("h_diff_l", "h_diff_l", 90, -60.0, 120);

    TH1F* h_ovlp = new TH1F("h_ovlp", "h_ovlp", 100, 0.0, 1.0);

    TH1F* h_vtx_x = new TH1F("h_vtx_x", "h_vtx_x", 50, -100 + tpc_x, 100 + tpc_x);
    TH1F* h_vtx_y = new TH1F("h_vtx_y", "h_vtx_y", 50, -100 + tpc_y, 100 + tpc_y);
    TH1F* h_vtx_z = new TH1F("h_vtx_z", "h_vtx_z", 50, -100 + tpc_z, 100 + tpc_z);


    h_num_ixn->SetTitle("Number of ML Reco Interactions Per Event;ML Reco Interactions;Count [Events/Interaction]");

    h_part_T->SetTitle("ML Reco Primary Track Energy;Energy [GeV];Count [Tracks/0.01 GeV]");
    h_part_p->SetTitle("ML Reco Primary Track Momentum;Momentum [GeV/cm];Count [Tracks/0.02 GeV/cm]");
    h_part_a->SetTitle("ML Reco Cosine of Primary Track Angle;Cosine of Track Angle with respect to Beam Direction;Count [Tracks/0.04]");
    h_part_a_rot->SetTitle("ML Reco Cosine of Primary Track Rotational Angle Projected onto Anode;Cosine of Track Rotational Angle Projected onto Anode;Count [Tracks/0.04]");
    h_part_a_incl->SetTitle("ML Reco Cosine of Primary Track Inclination Angle with respect to Anode;Cosine of Track Inclination Angle with respect to Anode;Count [Tracks/0.04]");
    h_part_l->SetTitle("ML Reco Primary Track Length;Length [cm];Count [Tracks/2 cm]");
    h_part_pdg->SetTitle("ML Reco Primary Track PDG;PDG;Count [Tracks/1]");

    h_true_T->SetTitle("True Primary Track Energy;Energy [GeV];Count [Tracks/0.01 GeV]");
    h_true_p->SetTitle("True Primary Track Momentum;Momentum [GeV/cm];Count [Tracks/0.02 GeV/cm]");
    h_true_a->SetTitle("True Cosine of Primary Track Angle;Cosine of Track Angle with respect to Beam Direction [#circ];Count [Tracks/0.04]");
    h_true_a_rot->SetTitle("True Cosine of Primary Track Rotational Angle Projected onto Anode;Cosine of Track Rotational Angle Projected onto Anode;Count [Tracks/0.04]");
    h_true_a_incl->SetTitle("True Cosine of Primary Track Inclination Angle with respect to Anode;Cosine of Track Inclination Angle with respect to Anode;Count [Tracks/0.04]");
    h_true_l->SetTitle("True Primary Track Length;Length [cm];Count [Tracks/2 cm]");
    h_true_pdg->SetTitle("True Primary Track PDG;PDG;Count [Tracks/1]");

    h_diff_T->SetTitle("Difference in ML Reco and True Primary Track Energy;Energy [GeV];Count [Tracks/0.01 GeV]");
    h_diff_p->SetTitle("Difference in ML Reco and True Primary Track Momentum;Momentum [GeV/cm];Count [Tracks/0.02 GeV/cm]");
    h_diff_a->SetTitle("Difference in ML Reco and True Cosine of Primary Track Angle;Cosine of Track Angle with respect to Beam Direction [#circ];Count [Tracks/0.04]");
    h_diff_a_rot->SetTitle("Difference in ML Reco and True Cosine of Primary Track Rotational Angle Projected onto Anode;Cosine of Track Rotational Angle Projected onto Anode;Count [Tracks/0.04]");
    h_diff_a_incl->SetTitle("Difference in ML Reco and True Cosine of Primary Track Inclination Angle with respect to Anode;Cosine of Track Inclination Angle with respect to Anode;Count [Tracks/0.04]");
    h_diff_l->SetTitle("Difference in ML Reco and True Primary Track Length;Length [cm];Count [Tracks/2 cm]");
    h_ovlp->SetTitle("Overlap of ML Reco and True Primary Track;Overlap;Count [Tracks/0.01]");

    h_vtx_x->SetTitle("Interaction ML Reco Vertex x-Position;x Position [cm];Count [Events/4 cm]");
    h_vtx_y->SetTitle("Interaction ML Reco Vertex y-Position;y Position [cm];Count [Events/4 cm]");
    h_vtx_z->SetTitle("Interaction ML Reco Vertex z-Position;z Position [cm];Count [Events/4 cm]");

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

        const auto num_ixn = sr->common.ixn.ndlp;
        h_num_ixn->Fill(num_ixn);

        for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
        {
            const auto& vtx = sr->common.ixn.dlp[ixn].vtx;

            //Get the truth interaction(s) corresponding to this reco interaction
            const auto& vec_truth_ixn = sr->common.ixn.dlp[ixn].truth;
            const auto& vec_overlap_ixn = sr->common.ixn.dlp[ixn].truthOverlap;

            if(vec_overlap_ixn.empty())
                continue;

            //Find the truth interaction with the largest overlap
            //auto result = std::max_element(vec_overlap_ixn.begin(), vec_overlap_ixn.end());
            //auto max_overlap = std::distance(vec_overlap_ixn.begin(), result);

            //Do this manually since SRProxy wraps all types and loses some functionality in the process
            //Maybe there is a way to get the object inside the Proxy, but this works for now
            double current_max = 0;
            unsigned int max_overlap = 0;
            for(unsigned int i = 0; i < vec_overlap_ixn.size(); ++i)
            {
                auto val = vec_overlap_ixn.at(i);
                if(val > current_max)
                {
                    current_max = val;
                    max_overlap = i;
                }
            }
            const auto truth_idx = vec_truth_ixn.at(max_overlap);
            const auto& truth_ixn = sr->mc.nu[truth_idx];

            //Put cuts on true interaction quantities here
            //For example, reject interactions not on argon
            //if(truth_ixn.targetPDG != 1000180400)
            //    continue;

            if(std::isfinite(vtx.x) && std::isfinite(vtx.y) && std::isfinite(vtx.z))
            {
                h_vtx_x->Fill(vtx.x);
                h_vtx_y->Fill(vtx.y);
                h_vtx_z->Fill(vtx.z);
            }

            auto reco_ixn_trks = 0;
            auto true_ixn_trks = 0;
            // Loop over particles in matched true interaction
            for(unsigned long ipart = 0; ipart < sr->mc.nu[truth_idx].prim.size(); ++ipart)
            {
                //Store current true particle for easier access
                const auto& true_part = sr->mc.nu[truth_idx].prim[ipart];
                if((true_part.pdg != 2212 and true_part.pdg != 321 and true_part.pdg != 211 and true_part.pdg != 13))
                    continue;
                ++true_ixn_trks;
            }
            // Loop over particles in the interaction
            for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
            {
                //Store current reco particle for easier access
                const auto& part = sr->common.ixn.dlp[ixn].part.dlp[ipart];

                //Put reco particle cuts here
                if((part.pdg != 2212 and part.pdg != 321 and part.pdg != 211 and part.pdg != 13) or not part.primary)
                    continue;

                //Put reco particle cuts here
                if((abs(abs(part.start.z)-64.538) < 1.0) and (abs(abs(part.end.z)-64.538)) < 1.0){
                    std::cout<<"Particle Start:"<<part.start.z<<" End:"<<part.end.z<<std::endl;
                    continue;
                }

                //Put reco particle cuts here
                if((abs(part.start.z+64.538) < 1.0) or (abs(part.end.z+64.538) < 1.0)){
                    std::cout<<"Particle Start:"<<part.start.z<<" End:"<<part.end.z<<std::endl;
                    continue;
                }

                
                //Get truth match(es) for this reco particle
                //Variables need to be wrapped in the Proxy object (sometimes)
                caf::Proxy<caf::SRTrueParticle>* truth_match = nullptr;
                const auto& vec_truth_id = part.truth;
                const auto& vec_overlap  = part.truthOverlap;

                //If the truth overlap vector is empty, then assume no truth match and skip
                if(vec_overlap.empty())
                {
                    //std::cout << "No truth match... skipping reco particle..." << std::endl;
                    continue;
                }


                //Find the truth particle with the largest overlap
                //auto result = std::max_element(vec_overlap.begin(), vec_overlap.end());
                //auto max_overlap = std::distance(vec_overlap.begin(), result);
                //auto truth_id = vec_truth_id.at(max_overlap);
                double current_max = 0;
                unsigned int max_overlap = 0;
                for(unsigned int i = 0; i < vec_overlap.size(); ++i)
                {
                    auto val = vec_overlap.at(i);
                    if(val > current_max)
                    {
                        current_max = val;
                        max_overlap = i;
                    }
                }

                //if(current_max < 0.25)
                //    continue;
                const auto& truth_id = vec_truth_id.at(max_overlap);

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

                ++reco_ixn_trks;

                //Finally get or calculate various reco/truth quantities
                auto pvec = TVector3(part.p.x, part.p.y, part.p.z);
                auto dir = TVector3(part.end.x, part.end.y, part.end.z) - TVector3(part.start.x, part.start.y, part.start.z);
                auto cos_angle = TMath::Cos(dir.Angle(beam_dir)); //calculate cosine of angle w.r.t neutrino beam direction
                dir.RotateY(-TMath::Pi()/2);
                auto cos_rot_anode_angle = TMath::Cos(dir.Theta()); //calculate cosine of track rotational angle (projection on anode)
                auto cos_incl_anode_angle = TMath::Cos(dir.Phi()); //calculate cosine of track inclination angle (off of anode)
                dir.RotateY(TMath::Pi()/2);
                //auto angle = pvec.Angle(beam_dir) * 180.0 / TMath::Pi(); //different method for the angle
                auto length = dir.Mag();

                auto true_pvec = TVector3(truth_match->p.px, truth_match->p.py, truth_match->p.pz);
                auto true_dir = TVector3(truth_match->end_pos.x, truth_match->end_pos.y, truth_match->end_pos.z)
                                - TVector3(truth_match->start_pos.x, truth_match->start_pos.y, truth_match->start_pos.z);
                auto true_cos_angle = TMath::Cos(true_dir.Angle(beam_dir));
                true_dir.RotateY(-TMath::Pi()/2);
                auto true_cos_rot_anode_angle = TMath::Cos(true_dir.Theta()); //calculate cosine of track rotational angle (projection on anode)
                auto true_cos_incl_anode_angle = TMath::Cos(true_dir.Phi()); //calculate cosine of track inclination angle (off of anode)
                true_dir.RotateY(TMath::Pi()/2);
                auto true_length_val = true_dir.Mag();

                auto T_diff = truth_match->p.E - part.E;
                auto p_diff = true_pvec.Mag() - pvec.Mag();
                auto length_diff = true_length_val - length;
                auto cos_angle_diff = true_cos_angle - cos_angle;

                dir.RotateY(-TMath::Pi()/2);
                true_dir.RotateY(-TMath::Pi()/2);

                // POPULATE: Record information in vectors (defined above) for tracks that
	            //           have passed all cuts
                
	            reco_energy.push_back(part.E);
                reco_p_x.push_back(part.p.x); 
                reco_p_y.push_back(part.p.y); 
                reco_p_z.push_back(part.p.z);
                reco_p_mag.push_back(pvec.Mag());
                reco_length.push_back(length);
                dir.RotateY(TMath::Pi()/2);
                reco_angle.push_back(dir.Angle(beam_dir));
                reco_angle_x.push_back(dir.Angle(x_plus_dir));
                reco_angle_y.push_back(dir.Angle(y_plus_dir));
                reco_angle_z.push_back(dir.Angle(z_plus_dir));
                dir.RotateY(-TMath::Pi()/2);
                reco_angle_rot.push_back(dir.Theta());
                reco_angle_incl.push_back(dir.Phi());
                reco_track_start_x.push_back(part.start.x);
                reco_track_start_y.push_back(part.start.y);
                reco_track_start_z.push_back(part.start.z);
                reco_track_end_x.push_back(part.end.x);
                reco_track_end_y.push_back(part.end.y);
                reco_track_end_z.push_back(part.end.z);
                reco_pdg.push_back(part.pdg);
                true_energy.push_back(truth_match->p.E);
                true_p_x.push_back(truth_match->p.px); 
                true_p_y.push_back(truth_match->p.py); 
                true_p_z.push_back(truth_match->p.pz);
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
                true_track_start_x.push_back(truth_match->start_pos.x);
                true_track_start_y.push_back(truth_match->start_pos.y);
                true_track_start_z.push_back(truth_match->start_pos.z);
                true_track_end_x.push_back(truth_match->end_pos.x);
                true_track_end_y.push_back(truth_match->end_pos.y);
                true_track_end_z.push_back(truth_match->end_pos.z);
                true_pdg.push_back(truth_match->pdg);
                true_ixn_charged_track_mult.push_back(true_ixn_trks);
                overlap.push_back(current_max);
                true_ixn_index.push_back(truth_idx);
                reco_ixn_index.push_back(ixn);
                spill_index.push_back(spill_num);
                file_index.push_back(file_num);


                //Reco kinematics cuts
                if(length > 60)
                    continue;

                h_part_T->Fill(part.E);
                h_part_p->Fill(pvec.Mag());
                h_part_a->Fill(cos_angle);
                h_part_a_rot->Fill(cos_rot_anode_angle);
                h_part_a_incl->Fill(cos_incl_anode_angle);
                h_part_l->Fill(length);
                h_part_pdg->Fill(part.pdg);

                h_true_T->Fill(truth_match->p.E);
                h_true_p->Fill(true_pvec.Mag());
                h_true_a->Fill(true_cos_angle);
                h_true_a_rot->Fill(true_cos_rot_anode_angle);
                h_true_a_incl->Fill(true_cos_incl_anode_angle);
                h_true_l->Fill(true_length_val);
                h_true_pdg->Fill(truth_match->pdg);

                h_diff_T->Fill(T_diff);
                h_diff_p->Fill(p_diff);
                h_diff_a->Fill(cos_angle_diff);
                h_diff_a_rot->Fill(true_cos_rot_anode_angle - cos_rot_anode_angle);
                h_diff_a_incl->Fill(true_cos_incl_anode_angle - cos_incl_anode_angle);
                h_diff_l->Fill(length_diff);

                h_ovlp->Fill(current_max);
            }
            // Loop over particles in the interaction again to load charged track multiplicity
            for(unsigned long ipart = 0; ipart < reco_ixn_trks; ++ipart)
            {
                reco_ixn_charged_track_mult.push_back(reco_ixn_trks);
            }
        }
    } //end of reco interactions

    } //end of file loop
    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    //Write plots to ROOT file, overwriting the ROOT file
    TFile* caf_output = new TFile("caf_plots.root", "recreate");
    h_num_ixn->Write();
    h_vtx_x->Write();
    h_vtx_y->Write();
    h_vtx_z->Write();
    h_part_T->Write();
    h_part_p->Write();
    h_part_a->Write();
    h_part_a_rot->Write();
    h_part_a_incl->Write();
    h_part_l->Write();
    h_part_pdg->Write();
    h_true_T->Write();
    h_true_p->Write();
    h_true_a->Write();
    h_true_a_rot->Write();
    h_true_a_incl->Write();
    h_true_l->Write();
    h_true_pdg->Write();
    h_diff_T->Write();
    h_diff_p->Write();
    h_diff_a->Write();
    h_diff_a_rot->Write();
    h_diff_a_incl->Write();
    h_diff_l->Write();
    h_ovlp->Write();
    caf_output->Close();

    //Write plots to PDF files, overwriting the previous files
    TCanvas* c_part_a_comp = new TCanvas("c_part_a_comp", "c_part_a_comp");
    c_part_a_comp->cd();
    h_part_a->SetTitle("ML Reco vs. Truth Cosine of Primary Track Angle;Cosine of Track Angle with respect to Beam Direction;Count [Tracks/0.04]");
    h_part_a->SetLineColor(kRed);
    h_part_a->Draw();
    h_true_a->SetTitle("ML Reco vs. Truth Cosine of Primary Track Angle;Cosine of Track Angle with respect to Beam Direction;Count [Tracks/0.04]");
    h_true_a->SetLineColor(kBlue);
    h_true_a->Draw("SAME");
    TLegend* leg_a = new TLegend(0.7, 0.7, 0.85, 0.85);
    leg_a->AddEntry(h_part_a, "ML Reco", "l");
    leg_a->AddEntry(h_true_a, "Truth", "l");
    leg_a->Draw();
    gStyle->SetOptStat(0); //turn off stats box
    c_part_a_comp->Update();
    c_part_a_comp->SaveAs("c_part_a_comp.pdf", "recreate");

    TCanvas* c_part_a_rot_comp = new TCanvas("c_part_a_rot_comp", "c_part_a_rot_comp");
    c_part_a_rot_comp->cd();
    h_part_a_rot->SetTitle("ML Reco vs. Truth Cosine of Primary Track Rotational Angle Projected onto Anode;Cosine of Track Rotational Angle Projected onto Anode;Count [Tracks/0.04]");
    h_part_a_rot->SetLineColor(kRed);
    h_part_a_rot->Draw();
    h_true_a_rot->SetTitle("ML Reco vs. Truth Cosine of Primary Track Rotational Angle Projected onto Anode;Cosine of Track Rotational Angle Projected onto Anode;Count [Tracks/0.04]");
    h_true_a_rot->SetLineColor(kBlue);
    h_true_a_rot->Draw("SAME");
    TLegend* leg_a_rot = new TLegend(0.7, 0.7, 0.85, 0.85);
    leg_a_rot->AddEntry(h_part_a_rot, "ML Reco", "l");
    leg_a_rot->AddEntry(h_true_a_rot, "Truth", "l");
    leg_a_rot->Draw();
    gStyle->SetOptStat(0); //turn off stats box
    c_part_a_rot_comp->Update();
    c_part_a_rot_comp->SaveAs("c_part_a_rot_comp.pdf", "recreate");

    TCanvas* c_part_a_incl_comp = new TCanvas("c_part_a_incl_comp", "c_part_a_incl_comp");
    c_part_a_incl_comp->cd();
    h_part_a_incl->SetTitle("ML Reco vs. Truth Cosine of Primary Track Inclination Angle with respect to Anode;Cosine of Track Inclination Angle with respect to Anode;Count [Tracks/0.04]");
    h_part_a_incl->SetLineColor(kRed);
    h_part_a_incl->Draw();
    h_true_a_incl->SetTitle("ML Reco vs. Truth Cosine of Primary Track Inclination Angle with respect to Anode;Cosine of Track Inclination Angle with respect to Anode;Count [Tracks/0.04]");
    h_true_a_incl->SetLineColor(kBlue);
    h_true_a_incl->Draw("SAME");
    TLegend* leg_a_incl = new TLegend(0.7, 0.7, 0.85, 0.85);
    leg_a_incl->AddEntry(h_part_a_incl, "ML Reco", "l");
    leg_a_incl->AddEntry(h_true_a_incl, "Truth", "l");
    leg_a_incl->Draw();
    gStyle->SetOptStat(0); //turn off stats box
    c_part_a_incl_comp->Update();
    c_part_a_incl_comp->SaveAs("c_part_a_incl_comp.pdf", "recreate");

    TCanvas* c_part_l_comp = new TCanvas("c_part_l_comp", "c_part_l_comp");
    c_part_l_comp->cd();
    h_part_l->SetTitle("ML Reco vs. Truth Primary Track Length;Length [cm];Count [Tracks/2 cm]");
    h_part_l->SetLineColor(kRed);
    h_part_l->Draw();
    h_true_l->SetTitle("ML Reco vs. Truth Primary Track Length;Length [cm];Count [Tracks/2 cm]");
    h_true_l->SetLineColor(kBlue);
    h_true_l->Draw("SAME");
    TLegend* leg_l = new TLegend(0.7, 0.7, 0.85, 0.85);
    leg_l->AddEntry(h_part_l, "ML Reco", "l");
    leg_l->AddEntry(h_true_l, "Truth", "l");
    leg_l->Draw();
    gStyle->SetOptStat(0); //turn off stats box
    c_part_l_comp->Update();
    c_part_l_comp->SaveAs("c_part_l_comp.pdf", "recreate");

    // Draw PDG Pie Chart for Truth Matches
    /*std::vector<char> true_pdg_codes;
    std::vector<double> true_pdg_counts;
    char counter; 
    for (int i = 1; i <= h_true_pdg->GetNbinsX(); i++) {
        if (h_true_pdg->GetBinContent(i) > 0) {
            true_pdg_codes.push_back(h_true_pdg->GetBin(i-1-10000));
            true_pdg_counts.push_back(h_true_pdg->GetBinContent(i));
            counter++;
            std::cout << h_true_pdg->GetBin(i-1-10000) << " " << h_true_pdg->GetBinContent(i) << std::endl;
        }
    }

    //const char** true_pdg_codes_labels = &true_pdg_codes;
    double* true_pdg_counts_arr = &true_pdg_counts[0];
    TCanvas* c_true_pdg_pie = new TCanvas("c_true_pdg_pie", "c_true_pdg_pie");
    c_true_pdg_pie->cd();
    TPie *p_true_pdg = new TPie("p_true_pdg", "Matched True Primary Track PDG", counter, true_pdg_counts_arr);
    p_true_pdg->SetCircle(.5,.5,.4);
    p_true_pdg->SetLabelsOffset(.01);
    //p_true_pdg->SetLabels(true_pdg_codes_labels);
    c_true_pdg_pie->Update();
    c_true_pdg_pie->SaveAs("c_true_pdg_pie.pdf", "recreate"); */

        // Output TTree file name
    std::string file_name = "track_reco_benchmark";

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
