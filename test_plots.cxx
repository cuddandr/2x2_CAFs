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

    TH1D* h_part_T = new TH1D("h_part_T", "h_part_T", 50, 0.0, 0.5);
    TH1D* h_part_p = new TH1D("h_part_p", "h_part_p", 50, 0.0, 1.0);

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
            for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
            {
                auto part = sr->common.ixn.dlp[ixn].part.dlp[ipart];
                //if(!part.contained)
                //    continue;

                if(part.pdg != 2212)
                    continue;

                float p = std::sqrt((part.p.x * part.p.x) + (part.p.y * part.p.y) + (part.p.z * part.p.z));

                h_part_T->Fill(part.E);
                h_part_p->Fill(p);
            }
        }
    }

    TFile* caf_output = new TFile("caf_plots.root", "recreate");
    h_part_T->Write();
    h_part_p->Write();
    caf_output->Close();

    std::cout << "Finished." << std::endl;
    return 0;
}
