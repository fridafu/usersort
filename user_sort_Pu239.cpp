#include "Event.h"
#include "Histogram1D.h"
#include "Histogram2D.h"
#include "IOPrintf.h"
#include "OfflineSorting.h"
#include "Parameters.h"
#include "ParticleRange.h"
#include "SiriusRoutine.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>


#define NDEBUG 1
#include "debug.h"

//! must be 0 to disable making DE:E matrices for each strip of all detectors
// Makes plots of type m_e_de_b0f0...b7f7
#define MAKE_INDIVIDUAL_E_DE_PLOTS 1

//! must be 0 to disable making time evolution plots
// This is to check, for example, if the elastic peak or a gamma peak 
// (which is more likely) moves during the experiment
#define MAKE_TIME_EVOLUTION_PLOTS 1

//! must be 0 to disable making CACTUS time-energy plots
// Makes plots of type m_nai_e_t_0...27
#define MAKE_CACTUS_TIME_ENERGY_PLOTS 1

//! must be 0 to disable checking particles
// If this is on, it checks that 
#define APPLY_PARTICLE_GATE 1

// ########################################################################

//! User sorting routine.
class UserXY : public SiriusRoutine {
public:
    UserXY();
    bool Sort(const Event& event);
    void CreateSpectra();
    bool Command(const std::string& cmd);
    
private:
    Histogram2Dp m_back, m_front, m_e_de_strip[8],/* m_e_de_dpEx0_strip[8],m_e_de_gamma_dpEx0_strip[8],*/m_e_de, m_e_de_thick, m_e_de_fiss, m_e_de_nofiss/*, m_e_de_gamma*/, m_e_de_fiss_nobg, m_e_de_fiss_bg;
#if defined(MAKE_INDIVIDUAL_E_DE_PLOTS) && (MAKE_INDIVIDUAL_E_DE_PLOTS>0)
    Histogram2Dp m_e_de_individual[8][8];
    Histogram1Dp h_ede_individual[8][8];
// add individual plot, zoomed in on d,p E_X=0 and for d,d*
//    Histogram1Dp h_e_dpEx0_individual[8][8];
// add individual plot, zoomed in on d,p E_X=0 and for d,d*
//    Histogram2Dp m_e_de_dpEx0_individual[8][8];
//    Histogram2Dp m_e_de_ddEx0_individual[8][8];
#endif /* MAKE_INDIVIDUAL_E_DE_PLOTS */
    Histogram2Dp m_nai_t, m_nai_e;
    Histogram2Dp m_alfna, m_alfna_bg;
    Histogram2Dp m_alfna_bg_nofiss, m_alfna_bg_fiss, m_alfna_fiss, m_alfna_nofiss;
    Histogram1Dp h_na_n, h_thick, /*h_thick_[8], h_thick_ind[8][8] ,*/h_ede, h_ede_r[8], h_ex_r[8], h_de_n, h_e_n;
    Histogram1Dp h_ex, h_ex_fiss;

#if defined(MAKE_CACTUS_TIME_ENERGY_PLOTS) && (MAKE_CACTUS_TIME_ENERGY_PLOTS>0)
    Histogram2Dp m_nai_e_t[28], m_nai_e_t_all, m_nai_e_t_c, m_siri_e_t[8], m_siri_e_t_all, m_ppac_e_t[4], m_ppac_e_t_all, m_ppac_e_t_c, m_nai_e_t_fiss[28], m_nai_e_t_all_fiss, m_nai_e_t_c_fiss, m_nai_e_t_fiss_bg[28], m_nai_e_t_all_fiss_bg, m_nai_e_t_c_fiss_bg;
#endif /* MAKE_CACTUS_TIME_ENERGY_PLOTS */

#if defined(MAKE_TIME_EVOLUTION_PLOTS) && (MAKE_TIME_EVOLUTION_PLOTS>0)
    Histogram2Dp m_nai_t_evol[28], m_nai_e_evol[28];
    Histogram2Dp m_e_evol[8], m_de_evol[8][8], m_ede_evol[8], m_ex_evol;
#endif /* MAKE_TIME_EVOLUTION_PLOTS */

private:
    //! Correction of CACTUS time for CACTUS energy.
    Parameter tnai_corr_enai;
    //! Correction of CACTUS time for SiRi back detector energy.
    Parameter tnai_corr_esi;

    //Gry trying to change the sort code so that we can read in separate parameters for correcting the PPAC time
    //This parameter is added to run_sunniva_deutrons.batch
    //! Correction of CACTUS time for SiRi back detector energy.
    Parameter tppac_corr_esi;

    //! Polynomials to calculate excitation energy from SiRi energy (E+DE).
    /*! Contains 8*3 coefficients. */
    Parameter ex_from_ede;
    //! Polynomials to make an empirical correction to the calculated excitation energy.
    /*! Contains 8*2 coefficients. */
    Parameter ex_corr_exp;
    //! Two rectangles to cut away SiRi noise/electrons.
    /*! Contains E-minimum 1, DE-minimum 1, E-minimum 2, DE-minimum 2. */
    Parameter ede_rect;
    //! Thickness centroid and second-order polynomial for SiRi E dependent thickness gate width.
    /*! Contains centroid, constant width, E-scaled width. */
    Parameter thick_range;
    //! The particle range data from zrange.
    ParticleRange particlerange;

    //! Apply energy corrections to CACTUS time.
    /*! \return corrected CACTUS time. */
    float tNaI(float t,    /*!< Uncorrected CACTUS time. */
               float Enai, /*!< Calibrated CACTUS energy in keV. */
               float Esi   /*!< Calibrated SiRi back energy in keV. */);

        //Declaration of function to correct PPAC time
    float tPpac(float t, /*!<Uncorrected CACTUS time.*/
                float Esi /*!Calibrated SiRi back energy in keV.*/);
  
    float range(float E /*!< particle energy in keV */)
        { return particlerange.GetRange( (int)E ); }
};

// ########################################################################

bool UserXY::Command(const std::string& cmd)
{
    std::istringstream icmd(cmd.c_str());

    std::string name;
    icmd >> name;
    
    if( name == "rangefile" ) {
        std::string filename;
        icmd >> filename;
        particlerange.Read( filename );
        return true;
    }
    return SiriusRoutine::Command(cmd);
}

// ########################################################################

UserXY::UserXY()
    : tnai_corr_enai ( GetParameters(), "tnai_corr_enai", 4 )
    , tnai_corr_esi  ( GetParameters(), "tnai_corr_esi", 4 )
    , tppac_corr_esi  ( GetParameters(), "tppac_corr_esi", 4 )
    , ex_from_ede    ( GetParameters(), "ex_from_ede", 8*3 )
    , ex_corr_exp    ( GetParameters(), "ex_corr_exp", 8*2 )
    , ede_rect       ( GetParameters(), "ede_rect", 4 )
    , thick_range    ( GetParameters(), "thick_range", 3 )
{
    ede_rect.Set( "500 250 30 500" );
    thick_range.Set( "130  13 0" );
}


// ########################################################################

void UserXY::CreateSpectra()
{
    const int max_e = 20000, max_de = 10000;
    m_back = Mat( "m_back", "back detector energies",
                  2000, 0, max_e, "E(Si) [keV]", 8, 0, 8, "detector nr." );
    m_front = Mat( "m_front", "front detector energies",
                   2000, 0, max_de, "#DeltaE(Si) [keV]", 64, 0, 64, "detector nr." );

#if defined(MAKE_INDIVIDUAL_E_DE_PLOTS) && (MAKE_INDIVIDUAL_E_DE_PLOTS>0)
    // try to change bin counts if not enough statistics for calibration
    for(int b=0; b<8; ++b ) {
        for(int f=0; f<8; ++f ) {
            m_e_de_individual[b][f] = 
                Mat( ioprintf("m_e_de_b%df%d", b, f), ioprintf("#DeltaE : E detector %d strip %d", b, f),
                     2000, 0, max_e, "E(Si) [keV]", 2000, 0, max_de, "#DeltaE(Si) [keV]" );
            h_ede_individual[b][f] =
                Spec( ioprintf("h_ede_b%df%d", b, f), ioprintf("E+#DeltaE detector %d strip %d", b, f),
                      2000, 0, max_e, "E+#DeltaE [keV]" );
// add infividual plot, zoomed in on d,p E_X=0 ... and for the d,d* reaction on Ex=0
//            h_e_dpEx0_individual[b][f] =
//                Spec( ioprintf("h_e_b%df%d", b, f), ioprintf("#Zoomed d,p DeltaE: E detector %d strip %d", b, f),
//                      200, 15000, 18500, "E [keV]" );
//            m_e_de_dpEx0_individual[b][f] = 
//                Mat( ioprintf("m_e_de_dpEx0_b%df%d", b, f), ioprintf("#Zoomed d,p DeltaE : E detector %d strip %d", b, f),
//                     200, 15000, 18500, "E(Si) [keV]", 200, 1000, 1800, "#DeltaE(Si) [keV]" );
//            m_e_de_ddEx0_individual[b][f] = 
//                Mat( ioprintf("m_e_de_ddEx0_b%df%d", b, f), ioprintf("#Zoomed d,d* DeltaE : E detector %d strip %d", b, f),
//                     200, 9600, 11000, "E(Si) [keV]", 200, 2800, 3600, "#DeltaE(Si) [keV]" );
//            h_thick_ind[b][f] = 
//                Spec( ioprintf("h_thick_b%df%d", b, f), ioprintf("apparent #DeltaE thickness : E detector %d strip %d", b, f),
//                     500, 0, 500, "#DeltaE 'thickness' [um]" );
        }
    }
#endif /* MAKE_INDIVIDUAL_E_DE_PLOTS */
    for(int f=0; f<8; ++f ) {
        m_e_de_strip[f] = Mat( ioprintf("m_e_de_f%d", f), ioprintf("E(NaI) : E(Si) strip %d", f),
                               2000, 0, max_e, "E(Si) [keV]", 2000, 0, max_de, "#DeltaE(Si) [keV]" );
    }
    // changed axis to compare the shift of the d,p' E_x=0 peak
//    for(int f=0; f<8; ++f ) {
//        m_e_de_dpEx0_strip[f] = Mat( ioprintf("m_e_de_dpEx0_f%d", f), ioprintf("E(NaI) : E(Si) strip %d", f),
//                               500, 15000, 18500, "E(Si) [keV]", 500, 1000, 1800, "#DeltaE(Si) [keV]" );
        // veto for gamma
//        m_e_de_gamma_dpEx0_strip[f] = Mat( ioprintf("m_e_de_gamma_dpEx0_f%d", f), ioprintf("E(NaI) : E(Si) strip %d", f),
//                               500, 15000, 18500, "E(Si) [keV]", 500, 1000, 1800, "#DeltaE(Si) [keV]" );
//    }

    m_e_de = Mat( "m_e_de", "#DeltaE : E for all detectors together",
                  500, 0, max_e, "E(Si) [keV]", 500, 0, max_de, "#DeltaE(Si) [keV]" );
    m_e_de_fiss = Mat( "m_e_de_fiss", "#DeltaE : E in coincidence with fission",
                 500, 0, max_e, "E(Si) [keV]", 500, 0, max_de, "#DeltaE(Si) [keV]" );
    m_e_de_nofiss = Mat( "m_e_de_nofiss", "#DeltaE : E veto for fission",
                 500, 0, max_e, "E(Si) [keV]", 500, 0, max_de, "#DeltaE(Si) [keV]" );
    m_e_de_thick = Mat( "m_e_de_thick", "#DeltaE : E for all detectors together, gated on thickness",
                        500, 0, max_e, "E(Si) [keV]", 500, 0, max_de, "#DeltaE(Si) [keV]" );
    m_e_de_fiss_bg = Mat( "m_e_de_fiss_bg", "#DeltaE : E in coincidence with fission, background",
                 500, 0, max_e, "E(Si) [keV]", 500, 0, max_de, "#DeltaE(Si) [keV]" );
    m_e_de_fiss_nobg = Mat( "m_e_de_fiss_nobg", "#DeltaE : E in coincidence with fission, w/o bg",
                 500, 0, max_e, "E(Si) [keV]", 500, 0, max_de, "#DeltaE(Si) [keV]" );
//    m_e_de_gamma = Mat( "m_e_de_gamma", "#DeltaE : E in coincidence with gamma",
//                 500, 0, max_e, "E(Si) [keV]", 500, 0, max_de, "#DeltaE(Si) [keV]" );
    
    m_nai_t = Mat( "m_nai_t", "t(NaI) matrix", 500, 0,  500, "? [a.u.]",     32,0,32, "det. id.");
    m_nai_e = Mat( "m_nai_e", "E(NaI) matrix", 2000, 0, 15000, "E(NaI) [keV]", 32,0,32, "det. id.");


    m_alfna = Mat( "m_alfna", "E(NaI) : E_{x}",
                   2000, -2000, 14000, "E(NaI) [keV]", 2000, -2000, 14000, "E_{x} [keV]" );
    m_alfna_fiss = Mat( "m_alfna_fiss", "E(NaI) : E_{x} in coincidence with fission",
                  2000, -2000, 14000, "E(NaI) [keV]", 2000, -2000, 14000, "E_{x} [keV]" );
    m_alfna_nofiss = Mat( "m_alfna_nofiss", "E(NaI) : E_{x} veto for fission",
                          2000, -2000, 14000, "E(NaI) [keV]", 2000, -2000, 14000, "E_{x} [keV]" );
    m_alfna_bg = Mat( "m_alfna_bg", "E(NaI) : E_{x} background",
                   2000, -2000, 14000, "E(NaI) [keV]", 2000, -2000, 14000, "E_{x} [keV]" );
    m_alfna_bg_fiss = Mat( "m_alfna_bg_fiss", "E(NaI) : E_{x} background with fission",
                  2000, -2000, 14000, "E(NaI) [keV]", 2000, -2000, 14000, "E_{x} [keV]" );
    m_alfna_bg_nofiss = Mat( "m_alfna_bg_nofiss", "E(NaI) : E_{x} background without fission",
                     2000, -2000, 14000, "E(NaI) [keV]", 2000, -2000, 14000, "E_{x} [keV]" );
    
    h_na_n = Spec("h_na_n", "NaI multiplicity", 32, 0, 32, "multiplicity");

    h_thick = Spec("h_thick", "apparent #DeltaE thickness", 500, 0, 500, "#DeltaE 'thickness' [um]");
    h_de_n = Spec("h_de_n", "#DeltaE multiplicity", 64, 0, 64, "multiplicity");
    h_e_n = Spec("h_e_n", "E multiplicity", 10, 0, 10, "multiplicity");
    

    for(int f=0; f<8; ++f ) {
        h_ede_r[f] = Spec(ioprintf("h_ede_f%d", f), ioprintf("E+#DeltaE ring %d", f),
                          2000, 0, max_e, "E+#DeltaE [keV]");
        //h_ede_r[f]->SetLineColor(f+1);

        h_ex_r[f] = Spec(ioprintf("h_ex_f%d", f), ioprintf("E_{x} ring %d", f),
                          2000, -2000, 14000, "E_{x} [keV]");
        //h_ex_r[f]->SetLineColor(f+1);
//        h_thick_[f] = Spec(ioprintf("h_thick_f%d", f), ioprintf("apparent #DeltaE thickness ring %d", f),
//                          500, 0, 500, "#DeltaE 'thickness' [um]");
    }
    h_ede = Spec("h_ede", "E+#DeltaE all detectors", 2000, 0, max_e, "E+#DeltaE [keV]");
    h_ex  = Spec("h_ex", "E_{x} all detectors", 2000, -2000, 14000, "E_{x} [keV]");
    h_ex_fiss  = Spec("h_ex_fiss", "E_{x} all detectors, in coincidence with fission", 2000, -2000, 14000, "E_{x} [keV]");

#if defined(MAKE_CACTUS_TIME_ENERGY_PLOTS) && (MAKE_CACTUS_TIME_ENERGY_PLOTS>0)
    const int max_enai = 12000;
    for(int n=0; n<28; ++n ) {
    m_nai_e_t[n] = Mat( ioprintf("m_nai_e_t_%02d", n), ioprintf("t : E NaI %d", n),
                            500, 0, max_enai, "E(NaI) [keV]", 500, 0, 500, "t(NaI) [a.u.]" );
    m_nai_e_t_fiss[n] = Mat( ioprintf("m_nai_e_t_fiss_%02d", n), ioprintf("t : E NaI %d", n, "coinc. PPAC"),
                            500, 0, max_enai, "E(NaI) [keV]", 500, 0, 500, "t(NaI) [a.u.]" );
    m_nai_e_t_fiss_bg[n] = Mat( ioprintf("m_nai_e_t_fiss_bg_%02d", n), ioprintf("t : E NaI %d", n, "coinc. PPAC bg."),
                            500, 0, max_enai, "E(NaI) [keV]", 500, 0, 500, "t(NaI) [a.u.]" );
    }
    m_nai_e_t_all = Mat( "m_nai_e_t", "t : E NaI all together",
                         500, 0, max_enai, "E(NaI) [keV]", 500, 0, 500, "t(NaI) [a.u.]" );
    m_nai_e_t_c   = Mat( "m_nai_e_t_c", "t : E NaI all together, corrected",
                         500, 0, max_enai, "E(NaI) [keV]", 500, 0, 500, "t(NaI) [a.u.]" );

    m_nai_e_t_all_fiss = Mat( "m_nai_e_t_fiss_all", "t : E NaI all together",
                         500, 0, max_enai, "E(NaI) [keV]", 500, 0, 500, "t(NaI) [a.u.], PPAC coincidence" );


    m_nai_e_t_c_fiss   = Mat( "m_nai_e_t_c_fiss", "t : E NaI all together, corr, PPAC coincidence",
                         500, 0, max_enai, "E(NaI) [keV]", 500, 0, 500, "t(NaI) [a.u.]" );

    m_nai_e_t_all_fiss_bg = Mat( "m_nai_e_t_fiss_bg_all", "t : E NaI all together, bg.",
                         500, 0, max_enai, "E(NaI) [keV]", 500, 0, 500, "t(NaI) [a.u.], PPAC coincidence" );

    m_nai_e_t_c_fiss_bg   = Mat( "m_nai_e_t_c_fiss_bg", "t : E NaI all together, corrected, PPAC coincidence bg.",
                         500, 0, max_enai, "E(NaI) [keV]", 500, 0, 500, "t(NaI) [a.u.]" );



    for(int n=0; n<8; ++n ) {
    m_siri_e_t[n] = Mat( ioprintf("m_siri_e_t_b%d", n), ioprintf("t(NaI) : E(Si) detector %d", n),
                             500, 0, max_e, "E(Si) [keV]", 500, 0, 500, "t(NaI) corr. [a.u.]" );
    }
    m_siri_e_t_all = Mat( "m_siri_e_t", "t(NaI) : E(Si) all detectors",
                          500, 0, max_e, "E(Si) [keV]", 500, 0, 500, "t(NaI) corr. [a.u.]" );

    

    for(int n=0; n<8; ++n ) {
    m_ppac_e_t[n] = Mat( ioprintf("m_ppac_e_t_b%d", n), ioprintf("t(PPAC) : E(Si) detector %d", n),
                         500, 0, max_e, "E(Si) [keV]", 500, 0, 500, "t(PPAC) [a.u.]" );
    }
    m_ppac_e_t_all = Mat( "m_ppac_e_t", "t(PPAC) : E(Si) all detectors",
                          500, 0, max_e, "E(Si) [keV]", 500, 0, 500, "t(PPAC) [a.u.]" );
    m_ppac_e_t_c = Mat( "m_ppac_e_t_c", "t(PPAC) : E(Si) all detectors, corrected",
                          500, 0, max_e, "E(Si) [keV]", 500, 0, 500, "t(PPAC) corr. [a.u.]" );



#endif /* MAKE_CACTUS_TIME_ENERGY_PLOTS */

#if defined(MAKE_TIME_EVOLUTION_PLOTS) && (MAKE_TIME_EVOLUTION_PLOTS>0)
    // time evolution plots
    const int MT = 4*24*3600;
    const int NT = 4*24;
    for(int n=0; n<28; ++n ) {
        m_nai_t_evol[n] = Mat( ioprintf("m_nai_t_evol_%02d", n), ioprintf("time : t NaI %d", n),
                            500, 0, 500, "t(NaI) [a.u.]", NT, 0, MT, "wall clock time [s]" );
        m_nai_e_evol[n] = Mat( ioprintf("m_nai_e_evol_%02d", n), ioprintf("time : e NaI %d", n),
                            500, -1000, max_enai-1000, "e(NaI) [keV]", NT, 0, MT, "wall clock time [s]" );
    }
    
    for(int b=0; b<8; ++b ) {
        m_e_evol[b] = Mat( ioprintf("m_e_evol_b%d", b), ioprintf("time : E detector %d", b),
                           500, 0, max_e, "E(Si) [keV]", NT, 0, MT, "wall clock time [s]" );
        for(int f=0; f<8; ++f ) {
            m_de_evol[b][f] = 
                Mat( ioprintf("m_de_evol_b%df%d", b, f), ioprintf("time : #DeltaE detector %d strip %d", b, f),
                     500, 0, max_de, "#DeltaE(Si) [keV]", NT, 0, MT, "wall clock time [s]" );
        }
        m_ede_evol[b] = Mat( ioprintf("m_ede_evol_f%d", b), ioprintf("time : E+#DeltaE ring %d", b),
                             500, 0, max_e, "E+#DeltaE(Si) [keV]", NT, 0, MT, "wall clock time [s]" );
    }
    m_ex_evol  = Mat("m_ex_evol", "time : E_{x} all detectors", 800, -2000, 14000, "E_{x} [keV]",
                     NT, 0, MT, "wall clock time [s]" );
#endif /* MAKE_TIME_EVOLUTION_PLOTS */
}

// ########################################################################

static float _rando = 0;
static float calib(unsigned int raw, float gain, float shift)
{
    return shift + (raw+_rando) * gain;
}

// ########################################################################

float UserXY::tNaI(float t, float Enai, float Esi)
{
    const float c = tnai_corr_enai[0] + tnai_corr_enai[1]/(Enai+tnai_corr_enai[2]) + tnai_corr_enai[3]*Enai;
    const float d = tnai_corr_esi [0] + tnai_corr_esi [1]/(Esi +tnai_corr_esi [2]) + tnai_corr_esi [3]*Esi;
    return t - c - d;
}

 // ########################################################################
//Function to correct the time of the PPAC
float UserXY::tPpac(float t, float Esi)
{  
    const float a = tppac_corr_esi [0] + tppac_corr_esi [1]/(Esi +tppac_corr_esi [2]) + tppac_corr_esi [3]*Esi;
    return t - a;
}
// ########################################################################


bool UserXY::Sort(const Event& event)
{
    _rando = drand48() - 0.5;
     
    
    // ..................................................
    // ALEXANDER's ORIGINAL ROUTINE
    // ..................................................
/*
    const unsigned int e_mini = (unsigned int)std::min(ede_rect[0], ede_rect[2]);
    unsigned int si_e_raw[8];
    for( int i=0; i<event.n_e; i++ ) {
        int id = event.e[i].chn;
        if( !(id&1) || id>= 16 )
            continue; // ignore guard rings

        id >>= 1; // detector number 0..7

        // only keep raw E here, we don't know which front detector
        // has fired, so we also don't know which coefficients to use
        // to calibrate the back detector
        const unsigned int raw = event.e[i].adc;
        if( raw >= e_mini )
            si_e_raw[id] = raw;
        else
            si_e_raw[id] = 0;

        // approximate calibration
        m_back->Fill( (int)calib( raw, gain_e[8*id], shift_e[8*id] ), id );
    }
    h_e_n->Fill(event.n_e);

    // ..................................................

    int si_goodcount = 0, dei=-1, ei=-1;
    float de = 0;

    // calibrate dE detectors, reject event if more than one over threshold
    for( int i=0; i<event.n_de; i++ ) {
        const int id   = event.de[i].chn;
        const int id_b = id / 8;
        const int id_f = id % 8;

        const unsigned int raw = event.de[i].adc;
        const float de_cal = calib( raw, gain_de[id], shift_de[id] );

        m_front->Fill( (int)de_cal, id );
        
        if( si_goodcount<2 && si_e_raw[id_b]>0 &&
            ( (si_e_raw[id_b]>ede_rect[0] && de_cal>ede_rect[1])
              || (si_e_raw[id_b]>ede_rect[2] && de_cal>ede_rect[3]) ) )
        {
            si_goodcount += 1;
            ei  = id_b;
            dei = id_f;
            de  = de_cal;
        }
    }
    
    h_de_n->Fill(event.n_de);
*/    
    
    // ......................................................//
    // MY MODIFIED SORTING ROUTINE
    // GIVES MUCH LESS "PILE-UP" TO THE LEFT OF ELASTIC PEAK
    // ......................................................//
    
    
    const int e_mini = 100;    // E threshold
    int si_e_raw[8];   // E energies
    int raw;
    
    int id_b;
    // ..................................................
    // E DETECTORS
    // Check E counters first, these give the master gate. 
    // Allow for only one E detector within the time range
    // (although there are many events with two and more!)
    // ..................................................
    
    for( int i=0; i<event.n_e; i++ ) {
    
    // ####### Fabio: If guard-rings fire, very bad!!!!!! ######## //
    // shame on "MY"//
//        if(event.n_e>1) // only one E with signal, else jump to next event

//            return true;
        
        int id = event.e[i].chn;
	
	//bad strip
	
        if( !(id&1) || id>= 16 )
            continue; // ignore guard rings. Guard rings in id 0,2,4,6,8,10,12,14, detectors in id 1,3,5,7,9,11,13,15
        
        id >>= 1; // detector id is now transformed into detector number 0..7
        id_b = id;  // keep the E ID
        
        // only keep raw E here, we don't know which front detector
        // has fired, so we also don't know which coefficients to use
        // to calibrate the back detector
        raw = event.e[i].adc;
        if( raw >= e_mini ) // check against E threshold
            si_e_raw[id] = raw; // if above, OK
        else{
            si_e_raw[id] = 0;   // if below, set energy to 0 and continue to next E counter
            continue;
        }
        
    }
    // approximate calibration
    m_back->Fill( (int)calib( raw, gain_e[8*id_b], shift_e[8*id_b] ), id_b );
    h_e_n->Fill(event.n_e);
    
    // ..................................................
    // DELTA E DETECTORS. 
    // NOW WE KNOW WHICH BACK DETECTOR HAS FIRED, SO WE
    // CAN CHECK THAT WE HAVE THE CORRECT FRONT DETECTOR.
    // CORRELATIONS:
    // id_back = 0, id_front = 0...7
    // id_back = 1, id_front = 8...15
    // id_back = 2, id_front = 16...23
    // id_back = 3, id_front = 24...31
    // id_back = 4, id_front = 32...39
    // id_back = 5, id_front = 40...47
    // id_back = 6, id_front = 48...55
    // id_back = 7, id_front = 56...63
    // ..................................................
    
    int si_goodcount = 0, dei=-1, ei=-1;
    float de = 0;
    int i_front_start = id_b*8;
    int i_front_stop  = (id_b*8) + 7;
    
    //    std::cout << " ............................................." << std::endl;
    //    std::cout << " Back detector fired: " << id_b <<  std::endl;
    //    std::cout << " Start front: " << i_front_start << ", stop front: " << i_front_stop << std::endl;
    //    std::cout << " Number of Delta E events:" << event.n_de << std::endl;
    for(int i=0;i<event.n_de;i++){
      
      
      
        const int id   = event.de[i].chn;
        //        std::cout << " Front ID fired: " << id << std::endl;
        if( id==29 ) //bad strip (id_b==3 && id_f==5) 
	    continue; 
	
        if(id<i_front_start || id>i_front_stop) // check id of front against the id of back - must match
            continue;
        
        //        std::cout << " Back ID: " << id_b << ", front strip " << id << ", energy front:" << raw << std::endl;
        
        const unsigned int raw = event.de[i].adc;
        const float de_cal = calib( raw, gain_de[id], shift_de[id] );
//        if(de_cal < 540)    // to exclude noise events, 106Cd exp.
        if(de_cal < 200)
            continue;
        
        m_front->Fill( (int)de_cal, id );   
        
        const int id_f = id % 8; // from here, the front strip has ID 0..7
        

        // Now ensure that only one strip has fired
        if(si_goodcount < 2){
            ei  = id_b;
            dei = id_f;
            de  = de_cal;
            si_goodcount += 1;
        }
    }
    
    
    h_de_n->Fill(event.n_de);  

    
    if( APPLY_PARTICLE_GATE && si_goodcount != 1 )
        // no detector above threshold, reject event
        return true;
 

    
    const float e  = calib( si_e_raw[ei], gain_e[8*ei+dei], shift_e[8*ei+dei] );
    const int e_int = int(e), de_int = int(de);

    // ..................................................

//****************************************************************************************************        
    // investigation for fission (SiRi)

    int fiss = 0;
    // int fissppac = 0;


    for( int j=0; j<event.n_na; j++ ) {
        
        const int ide = event.na[j].chn;
        
 //       if ( ide!=5 && ide!=13 && ide!=31 && ide!=32 )
        if ( ide!=4 && ide!=12 && ide!=30 && ide!=31 )
	  continue;
            
        
    //    const float na_e_f = calib( (int)event.na[j].adc, gain_na[idf], shift_na[idf] );
        
    //    const float na_t_f = calib( (int)event.na[j].tdc/8, gain_tna[idf], shift_tna[idf] );
        
        const float na_e_f = calib( (int)event.na[j].adc, gain_na[ide], shift_na[ide] );
        
        const float na_t_f = calib( (int)event.na[j].tdc/8, gain_tna[ide], shift_tna[ide] );      

//        if ( na_t_f>190 && na_t_f<220 && na_e_f>1195 && na_e_f<1225 ) fiss = 1;
// // Fabio: don't want energy requirement at the moment
//         if ( na_t_f>190 && na_t_f<220) fiss = 1;
//         if ( ppac_t_c>190 && ppac_t_c<220 &&  5e3 < e )   fissppac = 1; // select fission blob in tPPAC vs E_SiRi gate
//         if ( ppac_t_c>300 && ppac_t_c<330 &&  5e3 < e )   fissppac = 2; // added these to also see background fissions
    }
//****************************************************************************************************        
    
//****************************************************************************************************        

// gate on prompt d,d' peak
//    if (e_int<9800 ||  e_int>10250 ) 
//    return true;

    //E-dE matrix with veto for fission
    if( fiss==0 ) m_e_de_nofiss->Fill( e_int, de_int );
    
    //E-dE matrix only in case of fission
    if( fiss==1 ) m_e_de_fiss->Fill( e_int, de_int );
    
//****************************************************************************************************    
    
    
    // make DE:E matrices
#if defined(MAKE_INDIVIDUAL_E_DE_PLOTS) && (MAKE_INDIVIDUAL_E_DE_PLOTS>0) 
    m_e_de_individual[ei][dei]->Fill( e_int, de_int );
// for the zoomed version
//    m_e_de_dpEx0_individual[ei][dei]->Fill( e_int, de_int );
//    m_e_de_ddEx0_individual[ei][dei]->Fill( e_int, de_int );
#endif /* MAKE_INDIVIDUAL_E_DE_PLOTS */
    m_e_de_strip[dei]->Fill( e_int, de_int );
// make m_e_de_strip with certain axis, too,
//    m_e_de_dpEx0_strip[dei]->Fill( e_int, de_int );
    m_e_de->Fill( e_int, de_int );
    
    //const float thick = range(e+de)-range(e);
    const float thick = range(e+de)-range(e);
    h_thick->Fill( (int)thick );
//    h_thick_[dei]->Fill ( (int)thick );
//    h_thick_ind[ei][dei]-> Fill ( (int)thick );

    const float thick_dev = thick_range[1] + thick_range[2]*e;
    const bool have_pp = fabs(thick-thick_range[0])<thick_dev;
    if( APPLY_PARTICLE_GATE && !have_pp )
        return true;

    m_e_de_thick->Fill( e_int, de_int );
    const float ede = e+de;
    const int   ede_int = (int)ede;
    h_ede->Fill( ede_int );
    h_ede_r[dei]->Fill( ede_int );
#if defined(MAKE_INDIVIDUAL_E_DE_PLOTS) && (MAKE_INDIVIDUAL_E_DE_PLOTS>0)
    h_ede_individual[ei][dei]->Fill( ede_int );
//add individual h_e    
//    h_e_dpEx0_individual[ei][dei]->Fill( e_int );


#endif /* MAKE_INDIVIDUAL_E_DE_PLOTS */

    // fit of kinz Ex(E+DE)
    const float ex_theo = ex_from_ede[3*dei+0] + (ede)*(ex_from_ede[3*dei+1] + (ede)*ex_from_ede[3*dei+2]);
    //const float ex_theo = ex_from_ede.Poly(ede, 3*dei, 3);

    // make experimental corrections
    const float ex = ex_corr_exp[2*dei]+ex_corr_exp[2*dei+1]*ex_theo;
    const int   ex_int = (int)ex;

    h_ex->Fill( ex_int );
    h_ex_r[dei]->Fill( ex_int );

    if (fiss==1) {
    h_ex_fiss->Fill( ex_int );
    }
    // ..................................................

#if defined(MAKE_TIME_EVOLUTION_PLOTS) && (MAKE_TIME_EVOLUTION_PLOTS>0)
    const int timediff = Timediff(event);
    m_ex_evol->Fill( ex_int, timediff );
#endif /* MAKE_TIME_EVOLUTION_PLOTS */

    // ..................................................

    h_na_n->Fill(event.n_na);

    for( int i=0; i<event.n_na; i++ ) {
        const int id = event.na[i].chn;
	
  //      if( id == 15 || id == 9)
  //          continue;	

        
        //Gry har laget et ny "funksjon" UserXY::tPpac som leser inn NaI-tid og Siri-energi og bruker parameterene
        //tppac_corr_esi til å korrigere tiden.

        // const int   ppac_t_c = (int)tPpac(na_t,e);
    int fissppac = 0;


    // for( int j=0; j<event.n_na; j++ ) {
        // const int ide = event.na[j].chn;
        const int ide = event.na[i].chn;

 //       if ( ide!=5 && ide!=13 && ide!=31 && ide!=32 )
                int   ppac_t_c;
        if ( ide==4 || ide==12 || ide==30 || ide==31 ){

        const float na_e_f = calib( (int)event.na[i].adc, gain_na[ide], shift_na[ide] );
        
        const float na_t_f = calib( (int)event.na[i].tdc/8, gain_tna[ide], shift_tna[ide] );  

        const int   ppac_t_c = (int)tPpac(na_t_f,e); 
//        if ( na_t_f>190 && na_t_f<220 && na_e_f>1195 && na_e_f<1225 ) fiss = 1;
// Fabio: don't want energy requirement at the moment
        if ( na_t_f>190 && na_t_f<220) fiss = 1;
        if ( ppac_t_c>190 && ppac_t_c<230 &&  5e3 < e )   fissppac = 1; // select fission blob in tPPAC vs E_SiRi gate
        if ( ppac_t_c>300 && ppac_t_c<360 &&  5e3 < e )   fissppac = 2; // added these to also see background fissions
    // }
           }


        if( event.na[i].adc <= 0 )
            continue;

        if( event.na[i].tdc <= 0 )
            continue;
        const float na_t = calib( (int)event.na[i].tdc/8, gain_tna[id], shift_tna[id] ); 
        const int   na_t_int = (int)na_t;
        const float na_e = calib( (int)event.na[i].adc, gain_na[id], shift_na[id] );
        
        const int   na_t_c = (int)tNaI(na_t, na_e, e);

  //      std::cout << id << std::endl;
        

        
        const int   na_e_int = (int)na_e;
	
//	if (na_e_int<3500 || na_e_int>4000)
//	  continue;
	
        // PUT GATE ON PARTICLE PEAK 
	
/*	if (ex_int<5300 || ex_int>8000)
	  continue;
	
	if (ex_int<7400 && ex_int>5900)
	  continue;
*/

//	if (ex_int>6000)
//	  continue;
//        if(ex_int>1650 && ex_int<1918)    // 1st Ex in 28Si
        //if(ex_int>8705 && ex_int<9002)      // peak at 8904 keV in 28Si
        //if(ex_int>3550 && ex_int<3700)      // Gate in 57Fe(p,d)56Fe
//if (ex_int<10090 || ex_int>10575) //Gate on 9.1 MeV in Pu239, possibly the O16 1st ex?
//	  return true;
        m_nai_e->Fill( na_e_int, id );

        
        m_nai_t->Fill( na_t_int, id );

// //****************************************************************************************************        
        // // investigation for fission (CACTUS)
        // int fiss = 0;
        // for( int j=0; j<event.n_na; j++ ) {
            // const int ide = event.na[j].chn;
// //            if  ( ide!=5 && ide!=13 && ide!=31 && ide!=32 )
            // if  ( ide!=4 && ide!=12 && ide!=30 && ide!=31 )
                // continue;
            
// //            const float na_e_f = calib( (int)event.na[j].adc, gain_na[idf], shift_na[idf] );
            
// //            const float na_t_f = calib( (int)event.na[j].tdc/8, gain_tna[idf], shift_tna[idf] );
 // // Changed variable from idf -> ide
            
            
            // const float na_e_f = calib( (int)event.na[j].adc, gain_na[ide], shift_na[ide] );
            
            // const float na_t_f = calib( (int)event.na[j].tdc/8, gain_tna[ide], shift_tna[ide] );
            
            
// //            if ( na_t_f>190 && na_t_f<220 && na_e_f>1195 && na_e_f<1225 ) fiss = 1;
// // Fabio: don't want energy requirement at the moment
            // if ( na_t_f>190 && na_t_f<220) fiss = 1;
        // }
// //****************************************************************************************************     
    // Investigate once again 
// 
// //****************************************************************************************************        

	
	
#if defined(MAKE_CACTUS_TIME_ENERGY_PLOTS) && (MAKE_CACTUS_TIME_ENERGY_PLOTS>0)
//        if ( id!=5 && id!=13 && id!=31 && id!=32 && fiss==0) {   
        if ( id!=4 && id!=12 && id!=30 && id!=31 && fiss==0) {   
	    m_nai_e_t[id] ->Fill( na_e_int,  na_t_int );
            m_nai_e_t_all ->Fill( na_e_int,  na_t_int );
            m_nai_e_t_c   ->Fill( na_e_int,  na_t_c );
	    m_siri_e_t[ei]->Fill( e_int, na_t_c );
            m_siri_e_t_all->Fill( e_int, na_t_c );
        }
        
        if ( id!=4 && id!=12 && id!=30 && id!=31 && (ppac_t_c>190 && ppac_t_c<230 &&  5e3 < e)) {   
        m_nai_e_t_fiss[id] ->Fill( na_e_int,  na_t_int );
        m_nai_e_t_all_fiss ->Fill( na_e_int,  na_t_int );
        m_nai_e_t_c_fiss   ->Fill( na_e_int,  na_t_c );
        }

        if ( id!=4 && id!=12 && id!=30 && id!=31 && fissppac==2) {   
        m_nai_e_t_fiss_bg[id] ->Fill( na_e_int,  na_t_int );
        m_nai_e_t_all_fiss_bg ->Fill( na_e_int,  na_t_int );
        m_nai_e_t_c_fiss_bg   ->Fill( na_e_int,  na_t_c );
        }



        // if ( (id==4 || id==12 || id==30 || id==31) && fiss==1) {  //(do for any PPAC) AND (only when it is a fission)  
        // m_ppac_e_t[ei]->Fill( e_int, na_t_int );     // ppac are feeded in as a NaI signal, therefore we
        //     m_ppac_e_t_all->Fill( e_int, na_t_int ); // can use na_t_int as ppac times
        // }


        if ( (id==4 || id==12 || id==30 || id==31) ) {  //(do for any PPAC)
        m_ppac_e_t[ei]->Fill( e_int, na_t_int );     // ppac are feeded in as a NaI signal, therefore we
            m_ppac_e_t_all->Fill( e_int, na_t_int ); // can use na_t_int as ppac times
            m_ppac_e_t_c->Fill( e_int, ppac_t_c );
        }
       
#endif /* MAKE_CACTUS_TIME_ENERGY_PLOTS */

        // ..................................................
       
		/*** HERE COMES THE MAIN MATRIX FOR NaI ***/
		int weight = 1;

        //Particle-gamma matrix all together
//        if( id!=5 && id!=13 && id!=31 && id!=32 && na_t_c>185 && na_t_c<215 ) {
        if( id!=4 && id!=12 && id!=30 && id!=31 && na_t_c>190 && na_t_c<220) {
            m_alfna->Fill( na_e_int, ex_int, 1);
//        } else if(id!=5 && id!=13 && id!=31 && id!=32 && na_t_c>233 && na_t_c<263 ) {
        } else if(id!=4 && id!=12 && id!=30 && id!=31 && na_t_c>300 && na_t_c<330 ) {
            m_alfna->Fill( na_e_int, ex_int, -1);

            m_alfna_bg->Fill( na_e_int, ex_int );
            weight = -1;   
//Fabio: Here I also want a gamma to be present in my m_e_de_matrix''''
//            m_e_de_gamma_dpEx0_strip[dei]->Fill( e_int, de_int );
        }
        
//***************************************************************************************************        
        
        //Particle-gamma matrix with veto for fission
//        if( id!=5 && id!=13 && id!=31 && id!=32 && fiss==0 && na_t_c>201 && na_t_c<237 ) {
        if( id!=4 && id!=12 && id!=30 && id!=31 && fiss==0 && na_t_c>190 && na_t_c<220 ) {
                m_alfna_nofiss->Fill( na_e_int, ex_int, 1);
//            } else if( id!=5 && id!=13 && id!=31 && id!=32 && fiss==0 && na_t_c>111 && na_t_c<147 ) {
            } else if( id!=4 && id!=12 && id!=30 && id!=31 && fiss==0 && na_t_c>300&& na_t_c<330) {
                m_alfna_nofiss->Fill( na_e_int, ex_int, -1);
                m_alfna_bg_nofiss->Fill( na_e_int, ex_int );
            }

        //Particle-gamma matrix only in case of fission
//        if( id!=5 && id!=13 && id!=31 && id!=32 && fiss==1 && na_t_c>201 && na_t_c<237 ) {
        if( id!=4 && id!=12 && id!=30 && id!=31 && fiss==1 && na_t_c>190 && na_t_c<220 ) {
            m_alfna_fiss->Fill( na_e_int, ex_int, 1);
//        } else if( id!=5 && id!=13 && id!=31 && id!=32  && fiss==1 && na_t_c>111 && na_t_c<147 ) {
        } else if( id!=4 && id!=12 && id!=30 && id!=31  && fiss==1 && na_t_c>300 && na_t_c<330 ) {
            m_alfna_fiss->Fill( na_e_int, ex_int, -1);
            m_alfna_bg_fiss->Fill( na_e_int, ex_int );
        }
//****************************************************************************************************
//Particle-gamma matrix only in case of fission
//        if( id!=5 && id!=13 && id!=31 && id!=32 && fiss==1 && na_t_c>201 && na_t_c<237 ) {
//         if( id!=4 && id!=12 && id!=30 && id!=31 && fiss==1 && na_t_c>190 && na_t_c<220 ) {
//             //m_alfna_fiss->Fill( na_e_int, ex_int, 1);
//             m_e_de_fiss_nobg->Fill( e_int, de_int );
// //        } else if( id!=5 && id!=13 && id!=31 && id!=32  && fiss==1 && na_t_c>111 && na_t_c<147 ) {
//         } else if( id!=4 && id!=12 && id!=30 && id!=31  && fiss==1 && na_t_c>300 && na_t_c<330 ) {
//             m_e_de_fiss_nobg->Fill( e_int, de_int, -1);
//             m_e_de_fiss_bg->Fill( e_int, de_int );
//         }




 

#if defined(MAKE_TIME_EVOLUTION_PLOTS) && (MAKE_TIME_EVOLUTION_PLOTS>0)
        m_nai_e_evol[id]->Fill( na_e_int, timediff, weight );
        m_nai_t_evol[id]->Fill( na_t_c,   timediff );
#endif /* MAKE_TIME_EVOLUTION_PLOTS */
    }
#if defined(MAKE_TIME_EVOLUTION_PLOTS) && (MAKE_TIME_EVOLUTION_PLOTS>0)
    m_e_evol  [ei]     ->Fill( e_int,   timediff );
    m_de_evol [ei][dei]->Fill( de_int,  timediff );
    m_ede_evol[dei]    ->Fill( ede_int, timediff );
#endif /* MAKE_TIME_EVOLUTION_PLOTS */

    return true;
}

// ########################################################################
// ########################################################################
// ########################################################################

int main(int argc, char* argv[])
{
    return OfflineSorting::Run(new UserXY(), argc, argv );
}
