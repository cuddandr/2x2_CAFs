unsigned int is_true_signal(const caf::SRTrueInteraction& mc_nu)
{
    const unsigned int argon_pdg = 1000180400;
    const float tpc_border       = 65;
    const float wall_cut         = 8.0;

    // unsigned int cut_level = 0;
    const auto true_ixn = mc_nu;
    const auto true_vtx = mc_nu.vtx;

    if(true_ixn.targetPDG != argon_pdg)
        return 0;

    if(true_ixn.iscc == false)
        return 1;

    if(std::abs(true_ixn.pdg) != 14)
        return 2;

    bool cut_vtx_x
        = std::abs(true_vtx.x) > wall_cut && std::abs(true_vtx.x) < (tpc_border - wall_cut);
    bool cut_vtx_y = std::abs(true_vtx.y) < (tpc_border - wall_cut);
    bool cut_vtx_z
        = std::abs(true_vtx.z) > wall_cut && std::abs(true_vtx.z) < (tpc_border - wall_cut);
    if(!(cut_vtx_x && cut_vtx_y && cut_vtx_z))
        return 3;

    unsigned int npi = true_ixn.npip + true_ixn.npim + true_ixn.npi0;
    if(npi > 0)
        return 4;

    unsigned int nka = 0;
    for(unsigned int p = 0; p < true_ixn.prim.size(); ++p)
    {
        auto part    = true_ixn.prim[p];
        auto abs_pdg = std::abs(part.pdg);
        if(abs_pdg == 321)
            nka++;
        if(abs_pdg == 130 || abs_pdg == 310 || abs_pdg == 311)
            nka++;
    }
    if(nka > 0)
        return 5;

    return 6;
}
/*
unsigned int has_minerva_match(const caf::SRRecoParticle& reco_part,
                               const caf::SRMINERvAInt& mnv_ixn)
{
    float dx_part   = reco_part->end.x - reco_part->start.x;
    float dy_part   = reco_part->end.y - reco_part->start.y;
    float dz_part   = reco_part->end.z - reco_part->start.z;
    float len_part  = std::sqrt(dx_part * dx_part + dy_part * dy_part + dz_part * dz_part);
    float dirx_part = dx_part / len_part;
    float diry_part = dy_part / len_part;
    float dirz_part = dz_part / len_part;

    bool mnv_match = false;
    bool mnv_punch = false;

    const float mnv_dot_cut    = 0.9975;
    const float mnv_dir_cut    = 0.06;
    const float mnv_extrap_cut = 10;

    for(unsigned int imnv = 0; imnv < mnv_ixn.size(); ++imnv)
    {
        for(unsigned int mtrk = 0; mtrk < mnv_ixn[imnv].ntracks; ++mtrk)
        {
            const auto& mnv_track = mnv_ixn[imnv].tracks[mtrk];
            h_mnv_trk_z->Fill(mnv_track.start.z, mnv_track.end.z);

            if(mnv_track.start.z > 160)
            {
                float dx_mnv  = mnv_track.end.x - mnv_track.start.x;
                float dy_mnv  = mnv_track.end.y - mnv_track.start.y;
                float dz_mnv  = mnv_track.end.z - mnv_track.start.z;
                float len_mnv = std::sqrt(dx_mnv * dx_mnv + dy_mnv * dy_mnv + dz_mnv * dz_mnv);

                if(len_mnv < 10)
                    continue;

                float dirx_mnv = dx_mnv / len_mnv;
                float diry_mnv = dy_mnv / len_mnv;
                float dirz_mnv = dz_mnv / len_mnv;

                float extrap_z = mnv_track.start.z - reco_part->end.z;
                float extrap_y
                    = diry_part / dirz_part * extrap_z + reco_part->end.y - mnv_track.start.y;
                float extrap_x
                    = dirx_part / dirz_part * extrap_z + reco_part->end.x - mnv_track.start.x;

                float mnv_dot_product
                    = dirx_mnv * dirx_part + diry_mnv * diry_part + dirz_mnv * dirz_part;
                if(mnv_dot_cut < mnv_dot_product
                   && std::abs(extrap_y) < mnv_extrap_cut
                   && std::abs(extrap_x) < mnv_extrap_cut
                   && std::abs(std::acos(dirx_mnv) - std::acos(dirx_part)) < mnv_dir_cut
                   && std::abs(std::acos(diry_mnv) - std::acos(diry_part)) < mnv_dir_cut)
                {
                    mnv_match = true;
                    if(mnv_track.end.z > 300)
                        mnv_punch = true;
                }
            }
        }
    }

    if(mnv_punch)
        return 2;
    else if(mnv_match)
        return 1;
    else
        return 0;
}
*/
