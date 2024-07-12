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

    TH1F* h_num_ixn = new TH1F("h_num_ixn", "h_num_ixn", 10, 0, 10);
    TH1F* h_num_trk = new TH1F("h_num_trk", "h_num_trk", 10, 0, 10);

    TH1F* h_track_E = new TH1F("h_track_E", "h_track_E", 50, 0.0, 0.5);
    TH1F* h_track_L = new TH1F("h_track_L", "h_track_L", 100, 0, 1000);

    TH1F* h_start_x = new TH1F("h_start_x", "h_start_x", 200, -1000, 1000);
    TH1F* h_start_y = new TH1F("h_start_y", "h_start_y", 200, -1000, 1000);
    TH1F* h_start_z = new TH1F("h_start_z", "h_start_z", 200, -1000, 1000);

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

        const auto num_ixn = sr->nd.minerva.nixn;
        h_num_ixn->Fill(num_ixn);

        for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
        {
            const auto num_trk = sr->nd.minerva.ixn[ixn].ntracks;
            h_num_trk->Fill(num_trk);

            for(unsigned long itrk = 0; itrk < num_trk; ++itrk)
            {
                auto trk = sr->nd.minerva.ixn[ixn].tracks[itrk];

                h_track_E->Fill(trk.E);
                h_track_L->Fill(trk.len_cm);
            }
        }
    }
    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    TFile* caf_output = new TFile("caf_plots.root", "recreate");
    h_num_ixn->Write();
    h_num_trk->Write();
    h_track_E->Write();
    h_track_L->Write();
    caf_output->Close();

    std::cout << "Time elapsed: " << t_elapsed.count() << std::endl;
    std::cout << "Finished." << std::endl;
    return 0;
}
