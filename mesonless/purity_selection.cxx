#include "/cvmfs/dune.opensciencegrid.org/products/dune/duneanaobj/v03_02_01/include/duneanaobj/StandardRecord/StandardRecord.h"

#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "cuts.cxx"

int purity_selection(const std::string& file_list, bool verbose = false)
{
    const auto t_load{std::chrono::steady_clock::now()};
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
    TFile* events_output = new TFile("events_mesonless.root", "recreate");
    auto tree_events = caf_chain->CloneTree(0);

    std::cout << "Processing " << root_list.size() << " files..." << std::endl;
    const unsigned int argon_pdg = 1000180400;

    const auto beam_dir    = TVector3(0.0, -0.05836, 1.0);
    const float tpc_border = 65;
    const float wall_cut   = 8.0;

    unsigned int total_signal_ixn = 0;
    unsigned int total_selected = 0;
    unsigned int signal_selected = 0;

    TH1I* h_num_true    = new TH1I("h_num_true", "h_num_true", 7, 0, 7);
    TH1I* h_tru_reac    = new TH1I("h_tru_reac", "h_tru_reac", 11, 0, 11);
    TH2F* h_true_muon   = new TH2F("h_true_muon", "h_true_muon", 80, 0, 8, 60, 0, 180);
    TH2F* h_true_proton = new TH2F("h_true_proton", "h_true_proton", 80, 0, 4, 60, 0, 180);

    TH1I* h_num_reco    = new TH1I("h_num_reco", "h_num_reco", 5, 0, 5);
    TH1I* h_sel_reac    = new TH1I("h_sel_reac", "h_sel_reac", 11, 0, 11);
    TH2F* h_reco_vtx    = new TH2F("h_reco_vtx", "h_reco_vtx", 60, -60, 60, 60, -60, 60);
    TH2F* h_reco_muon   = new TH2F("h_reco_muon", "h_reco_muon", 40, 0, 2, 60, 0, 180);
    TH2F* h_reco_proton = new TH2F("h_reco_proton", "h_reco_proton", 40, 0, 2, 60, 0, 180);

    TH2F* h_mnv_trk_z = new TH2F("h_mnv_trk_z", "h_mnv_trk_z", 180, -330, 330, 180, -330, 330);

    auto sr = new caf::StandardRecord;
    caf_chain->SetBranchAddress("rec", &sr);

    const unsigned long nspills = caf_chain->GetEntries();
    const unsigned int incr     = nspills / 10;
    const auto t_loop{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_intro{t_loop - t_load};

    std::cout << "Time loading: " << t_intro.count() << std::endl;
    std::cout << "Looping over " << nspills << " entries/spills..." << std::endl;

    const auto t_start{std::chrono::steady_clock::now()};
    for(unsigned long i = 0; i < nspills; ++i)
    {
        caf_chain->GetEntry(i);

        if(i % incr == 0)
            std::cout << "Spill #: " << i << std::endl;

        // Select all truth interactions
        for(unsigned long t = 0; t < sr->mc.nu.size(); ++t)
        {
            auto true_ixn = sr->mc.nu[t];
            auto true_vtx = sr->mc.nu[t].vtx;

            /*
            h_num_true->Fill(0.5);
            if(true_ixn.targetPDG != argon_pdg)
                continue;

            h_num_true->Fill(1.5);
            if(true_ixn.iscc == false)
                continue;

            h_num_true->Fill(2.5);
            if(std::abs(true_ixn.pdg) != 14)
                continue;

            h_num_true->Fill(3.5);
            bool cut_vtx_x
                = std::abs(true_vtx.x) > wall_cut && std::abs(true_vtx.x) < (tpc_border - wall_cut);
            bool cut_vtx_y = std::abs(true_vtx.y) < (tpc_border - wall_cut);
            bool cut_vtx_z
                = std::abs(true_vtx.z) > wall_cut && std::abs(true_vtx.z) < (tpc_border - wall_cut);
            if(!(cut_vtx_x && cut_vtx_y && cut_vtx_z))
                continue;

            h_num_true->Fill(4.5);
            unsigned int npi = true_ixn.npip + true_ixn.npim + true_ixn.npi0;
            if(npi > 0)
                continue;

            unsigned int nka = 0;
            for(unsigned int p = 0; p < true_ixn.prim.size(); ++p)
            {
                auto part    = true_ixn.prim[p];
                auto abs_pdg = std::abs(part.pdg);
                if(abs_pdg == 321)
                    nka++;
                if(abs_pdg == 130 || abs_pdg == 310)
                    nka++;
            }

            h_num_true->Fill(5.5);
            if(nka > 0)
                continue;

            h_num_true->Fill(6.5);
            */

            auto cut_level = is_true_signal(true_ixn);
            for(unsigned int c = 0; c <= cut_level; ++c)
                h_num_true->Fill(c + 0.5);

            if(cut_level < 6)
                continue;

            total_signal_ixn++;
            for(unsigned int k = 0; k < true_ixn.prim.size(); ++k)
            {
                auto part = true_ixn.prim[k];

                //auto E = part.p.E;
                auto p = part.p.Mag();

                auto dir   = TVector3(part.end_pos) - TVector3(part.start_pos);
                auto angle = part.p.Vect().Angle(beam_dir) * 180.0 / TMath::Pi();
                auto len   = dir.Mag();

                if(std::abs(part.pdg) == 13)
                    h_true_muon->Fill(p, angle);

                if(std::abs(part.pdg) == 2212)
                    h_true_proton->Fill(p, angle);
            }
        } // end truth interactions

        // Select reco interactions
        const auto num_ixn = sr->common.ixn.ndlp;
        for(unsigned long ixn = 0; ixn < num_ixn; ++ixn)
        {
            const auto& vtx = sr->common.ixn.dlp[ixn].vtx;
            h_num_reco->Fill(0.5);

            bool cut_vtx_x
                = std::abs(vtx.x) > wall_cut && std::abs(vtx.x) < (tpc_border - wall_cut);
            bool cut_vtx_y = std::abs(vtx.y) < (tpc_border - wall_cut);
            bool cut_vtx_z
                = std::abs(vtx.z) > wall_cut && std::abs(vtx.z) < (tpc_border - wall_cut);
            if(!(cut_vtx_x && cut_vtx_y && cut_vtx_z))
                continue;

            h_num_reco->Fill(1.5);
            const auto& vec_part   = sr->common.ixn.dlp[ixn].part.dlp;
            unsigned int num_muons = 0;
            unsigned int num_pions = 0;
            unsigned int num_kaons = 0;

            const caf::SRRecoParticle* reco_muon = nullptr;
            for(const auto& p : vec_part)
            {
                if(p.primary == false)
                    continue;

                if(std::abs(p.pdg) == 13)
                {
                    num_muons++;
                    reco_muon = &p;
                }

                if(std::abs(p.pdg) == 211)
                    num_pions++;

                if(std::abs(p.pdg) == 321)
                    num_kaons++;
            }

            if((num_muons != 1) or (num_pions > 0) or (num_kaons > 0))
                continue;

            h_num_reco->Fill(2.5);
            h_reco_vtx->Fill(vtx.x, vtx.z);

            float dx_muon   = reco_muon->end.x - reco_muon->start.x;
            float dy_muon   = reco_muon->end.y - reco_muon->start.y;
            float dz_muon   = reco_muon->end.z - reco_muon->start.z;
            float len_muon  = std::sqrt(dx_muon * dx_muon + dy_muon * dy_muon + dz_muon * dz_muon);
            float dirx_muon = dx_muon / len_muon;
            float diry_muon = dy_muon / len_muon;
            float dirz_muon = dz_muon / len_muon;

            bool mnv_match = false;
            bool mnv_punch = false;

            const float mnv_dot_cut    = 0.9975;
            const float mnv_dir_cut    = 0.06;
            const float mnv_extrap_cut = 10;
            for(unsigned int imnv = 0; imnv < sr->nd.minerva.ixn.size(); ++imnv)
            {
                for(unsigned int mtrk = 0; mtrk < sr->nd.minerva.ixn[imnv].ntracks; ++mtrk)
                {
                    const auto& mnv_track = sr->nd.minerva.ixn[imnv].tracks[mtrk];
                    h_mnv_trk_z->Fill(mnv_track.start.z, mnv_track.end.z);

                    if(mnv_track.start.z > 160)
                    {
                        float dx_mnv = mnv_track.end.x - mnv_track.start.x;
                        float dy_mnv = mnv_track.end.y - mnv_track.start.y;
                        float dz_mnv = mnv_track.end.z - mnv_track.start.z;
                        float len_mnv
                            = std::sqrt(dx_mnv * dx_mnv + dy_mnv * dy_mnv + dz_mnv * dz_mnv);

                        if(len_mnv < 10)
                            continue;

                        float dirx_mnv = dx_mnv / len_mnv;
                        float diry_mnv = dy_mnv / len_mnv;
                        float dirz_mnv = dz_mnv / len_mnv;

                        float extrap_z = mnv_track.start.z - reco_muon->end.z;
                        float extrap_y = diry_muon / dirz_muon * extrap_z + reco_muon->end.y
                                         - mnv_track.start.y;
                        float extrap_x = dirx_muon / dirz_muon * extrap_z + reco_muon->end.x
                                         - mnv_track.start.x;

                        float mnv_dot_product = dirx_mnv * dirx_muon + diry_mnv * diry_muon + dirz_mnv * dirz_muon;

                        if(mnv_dot_cut < mnv_dot_product
                           && std::abs(extrap_y) < mnv_extrap_cut
                           && std::abs(extrap_x) < mnv_extrap_cut
                           && std::abs(std::acos(dirx_mnv) - std::acos(dirx_muon)) < mnv_dir_cut
                           && std::abs(std::acos(diry_mnv) - std::acos(diry_muon)) < mnv_dir_cut)
                        {
                            mnv_match = true;
                            if(mnv_track.end.z > 300)
                                mnv_punch = true;
                        }
                    }
                }
            }

            if(mnv_match == false)
                continue;
            h_num_reco->Fill(3.5);

            if(mnv_punch == false)
                continue;
            h_num_reco->Fill(4.5);

            total_selected++;
            tree_events->Fill();

            // Select reconstructed particles
            for(unsigned long ipart = 0; ipart < sr->common.ixn.dlp[ixn].part.dlp.size(); ++ipart)
            {
                auto part    = sr->common.ixn.dlp[ixn].part.dlp[ipart];
                auto abs_pdg = std::abs(part.pdg);
                if(part.primary == false)
                    continue;

                // float p = std::sqrt((part.p.x * part.p.x) + (part.p.y * part.p.y) + (part.p.z *
                // part.p.z));
                float p = part.p.Mag();

                auto dir = TVector3(part.end) - TVector3(part.start);
                // auto angle = dir.Angle(beam_dir) * 180.0 / TMath::Pi();
                auto angle = TVector3(part.p).Angle(beam_dir) * 180.0 / TMath::Pi();
                auto len   = dir.Mag();

                if(abs_pdg == 13)
                {
                    h_reco_muon->Fill(p, angle);
                }

                if(abs_pdg == 2212)
                {
                    h_reco_proton->Fill(p, angle);
                }
            }

            const auto& vec_truth_ixn = sr->common.ixn.dlp[ixn].truth;
            const auto& vec_overlap = sr->common.ixn.dlp[ixn].truthOverlap;
            if(vec_overlap.empty())
                continue;

            auto result = std::max_element(vec_overlap.begin(), vec_overlap.end());
            auto max_overlap = std::distance(vec_overlap.begin(), result);
            auto truth_id = vec_truth_ixn.at(max_overlap);

            /* const auto idx = vec_truth_ixn.at(truth_id); */
            const auto& truth_ixn = sr->mc.nu[truth_id];
            h_sel_reac->Fill(truth_ixn.mode + 0.5);
            /*
            unsigned int npi = truth_ixn.npip + truth_ixn.npim + truth_ixn.npi0;
            if(npi > 0)
                continue;
            */
            auto cut_level = is_true_signal(truth_ixn);
            if(cut_level < 6)
                continue;
            signal_selected++;

            h_tru_reac->Fill(truth_ixn.mode + 0.5);

            if(verbose)
            {
                std::cout << "------------" << std::endl;
                std::cout << "File idx " << caf_chain->GetFile()->GetName() << std::endl;
                std::cout << "Spll idx " << i << std::endl;
                std::cout << "Reco idx " << ixn << std::endl;
                std::cout << "Meta evt " << sr->meta.nd_lar.event << std::endl;
                std::cout << "True idx " << truth_id << std::endl;
                std::cout << "VtxID : " << truth_ixn.id << " w/ overlap " << vec_overlap.at(max_overlap) << std::endl;
                std::cout << "Nu PDG: " << truth_ixn.pdg << std::endl;
                std::cout << "React : " << std::boolalpha << truth_ixn.iscc << " | " << truth_ixn.mode << std::endl;
                std::cout << "N prtn: " << truth_ixn.nproton << std::endl;
                std::cout << "N neut: " << truth_ixn.nneutron << std::endl;
                std::cout << "N pip : " << truth_ixn.npip << std::endl;
                std::cout << "N pim : " << truth_ixn.npim << std::endl;
                std::cout << "N pi0 : " << truth_ixn.npi0 << std::endl;

                std::cout << "Nprim : " << truth_ixn.nprim << std::endl;
                std::cout << "Nsec  : " << truth_ixn.nsec  << std::endl;
            }
        } // end reco interactions
    }
    const auto t_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> t_elapsed{t_end - t_start};

    events_output->Write();
    /* events_output->Close(); */

    TFile* caf_output = new TFile("purity_plots.root", "recreate");
    caf_output->cd();
    h_num_true->Write();
    h_tru_reac->Write();
    h_true_muon->Write();
    h_true_proton->Write();
    h_num_reco->Write();
    h_sel_reac->Write();
    h_reco_vtx->Write();
    h_reco_muon->Write();
    h_reco_proton->Write();
    h_mnv_trk_z->Write();

    caf_output->Close();
    float efficiency = (float)signal_selected / (float)total_signal_ixn;
    float purity = (float)signal_selected / (float)total_selected;
    std::cout << "Purity   : " << signal_selected << " / " << total_selected << " ~ " << purity << std::endl;
    std::cout << "Efficency: " << signal_selected << " / " << total_signal_ixn << " ~ " << efficiency << std::endl;
    std::cout << "Time elapsed: " << t_elapsed.count() << std::endl;
    std::cout << "Finished." << std::endl;
    return 0;
}

int main(int argc, char** argv)
{
    std::string filelist = argv[1];
    unsigned int stat = purity_selection(filelist, false);
    return stat;
}
