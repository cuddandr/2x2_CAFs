struct Inspector
{
    TChain* caf_chain;
    unsigned long curr_evt;
    unsigned long nspills;

    TVector3 beam_dir;

    caf::StandardRecord* sr;

    Inspector(const std::string& file_list);
    ~Inspector();

    void NextEvent();
    void LoadEvent(unsigned long evt);
    void ListInteractions();
    void ListTrueParticles();
    void ListRecoParticles();
};

Inspector::Inspector(const std::string& file_list)
    : caf_chain(nullptr), sr(nullptr), curr_evt(0), nspills(0)
{
    std::vector<std::string> root_list;
    std::ifstream fin(file_list, std::ios::in);
    if(!fin.is_open())
    {
        std::cerr << "Failed to open " << file_list << std::endl;
        std::cerr << "Exiting" << std::endl;
        //return 111;
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
        //return 121;
    }

    caf_chain = new TChain("cafTree");
    for(const auto& file : root_list)
    {
        std::cout << "Adding " << file << " to TChain." << std::endl;
        caf_chain->Add(file.c_str());
    }
    std::cout << "Processing " << root_list.size() << " files..." << std::endl;

    caf_chain->SetBranchAddress("rec", &sr);

    nspills = caf_chain->GetEntries();
    std::cout << "Total " << nspills << " spills." << std::endl;
    std::cout << "Loading entry 0 by default..." << std::endl;
    caf_chain->GetEntry(curr_evt);

    std::cout << "Initializing some constants..." << std::endl;
    beam_dir = TVector3(0.0, -0.05836, 1.0);

    std::cout << std::fixed << std::setprecision(4);
}

Inspector::~Inspector()
{
    if(sr != nullptr)
        delete sr;

    if(caf_chain != nullptr)
        delete caf_chain;
}

void Inspector::NextEvent()
{
    const auto evt = curr_evt + 1;
    if(evt < nspills)
    {
        curr_evt = evt;
        caf_chain->GetEntry(evt);
    }
    else
        std::cout << "Next event out of bounds." << std::endl;
}

void Inspector::LoadEvent(unsigned long evt)
{
    if(evt < nspills)
    {
        curr_evt = evt;
        caf_chain->GetEntry(evt);
    }
    else
        std::cout << evt << " out of bounds." << std::endl;
}

void Inspector::ListInteractions()
{
    std::cout << "Spill " << curr_evt << std::endl;
    const auto num_ixn = sr->common.ixn.ndlp;
    for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
    {
        std::cout << "-------------------" << std::endl;
        std::cout << "Reco interaction " << ixn << std::endl;
        const auto& vec_truth_ixn = sr->common.ixn.dlp[ixn].truth;
        const auto& vec_overlap = sr->common.ixn.dlp[ixn].truthOverlap;

        for(unsigned long t = 0; t < vec_truth_ixn.size(); ++t)
        {
            auto idx = vec_truth_ixn.at(t);
            const auto& truth_ixn = sr->mc.nu[idx];
            std::cout << "Truth Vtx ID " << truth_ixn.id << " with overlap " << vec_overlap.at(t) << std::endl;
            std::cout << "Target: " << truth_ixn.targetPDG << std::endl;
            std::cout << "isRock: " << std::boolalpha << (truth_ixn.id < 1E9 ? false : true) << std::endl;
            std::cout << "isCC  : " << std::boolalpha << truth_ixn.iscc << std::endl;
            std::cout << "React : " << truth_ixn.mode << std::endl;
            std::cout << "PDG   : " << truth_ixn.pdg << " | " << "E = " << truth_ixn.E << std::endl;
            std::cout << "N prtn: " << truth_ixn.nproton << std::endl;
            std::cout << "N neut: " << truth_ixn.nneutron << std::endl;
            std::cout << "N pip : " << truth_ixn.npip << std::endl;
            std::cout << "N pim : " << truth_ixn.npim << std::endl;
            std::cout << "N pi0 : " << truth_ixn.npi0 << std::endl;

            std::cout << "Nprim : " << truth_ixn.nprim << std::endl;
            std::cout << "Nsec  : " << truth_ixn.nsec  << std::endl;
            std::cout << std::endl;
        }
    }
    std::cout << "-------------------" << std::endl;
}

void Inspector::ListTrueParticles()
{
    const auto num_ixn = sr->common.ixn.ndlp;
    for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
    {
        std::cout << "-------------------" << std::endl;
        std::cout << "Reco interaction " << ixn << std::endl;
        const auto& vec_truth_ixn = sr->common.ixn.dlp[ixn].truth;
        const auto& vec_overlap = sr->common.ixn.dlp[ixn].truthOverlap;

        for(unsigned long t = 0; t < vec_truth_ixn.size(); ++t)
        {
            auto idx = vec_truth_ixn.at(t);
            const auto& truth_ixn = sr->mc.nu[idx];
            std::cout << "Truth Vtx ID " << truth_ixn.id << " with overlap " << vec_overlap.at(t) << std::endl;

            std::cout << "Truth Primaries: " << std::endl;
            std::cout << "\tG4ID | PDG | P_mag" << std::endl;
            for(unsigned int p = 0; p < truth_ixn.nprim; ++p)
            {
                auto q = truth_ixn.prim[p];
                std::cout << "\t" << q.G4ID << " | " << q.pdg << " | " << q.p.Mag() << std::endl;
                //std::cout << q.start_pos.X() << ", " << q.start_pos.Y() << ", " << q.start_pos.Z() << std::endl;
                //std::cout << q.start_pos.Mag() << std::endl;
                //std::cout << q.end_pos.Mag() << std::endl;
            }

            std::cout << "Truth Secondaries: " << std::endl;
            std::cout << "\tG4ID | PDG | P_mag" << std::endl;
            for(unsigned int s = 0; s < truth_ixn.nsec; ++s)
            {
                auto q = truth_ixn.sec[s];
                std::cout << "\t" << q.G4ID << " | " << q.pdg << " | " << q.p.Mag() << std::endl;
                //std::cout << q.start_pos.X() << ", " << q.start_pos.Y() << ", " << q.start_pos.Z() << std::endl;
                //std::cout << q.start_pos.Mag() << std::endl;
                //std::cout << q.end_pos.Mag() << std::endl;
            }
        }
    }
}

void Inspector::ListRecoParticles()
{
    const auto num_ixn = sr->common.ixn.ndlp;
    for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
    {
        const auto& vec_part = sr->common.ixn.dlp[ixn].part.dlp;

        std::cout << "Reco interaction " << ixn << std::endl;
        std::cout << "Reco particles: ";
        std::cout << "[\n";
        std::cout << "PDG | momentum | angle | length\n";
        for(const auto& part : vec_part)
        {
            if(part.primary == false)
                continue;

            float pmag = std::sqrt((part.p.x * part.p.x) + (part.p.y * part.p.y) + (part.p.z * part.p.z));
            float angle = TVector3(part.p).Angle(beam_dir) * 180.0 / TMath::Pi();
            auto dir = TVector3(part.end) - TVector3(part.start);
            auto len = dir.Mag();
            std::cout << std::setw(4) << part.pdg << " | "
                      << std::setw(9) << pmag << " | "
                      << std::setw(8) << angle << " | "
                      << std::setw(8) << len << std::endl;
        }
        std::cout << "]\n";
    }
}
