// Script to plot all Y data
// with theoretical models
// Cecilie, 20 Feb 2018
// NEW data from Gry 25 Jan 2019 - new fit from 4 March 2019

// Function for E1 strength 
double FitFunctionE1(double *x, double *par){
    const double Pi = 3.14159;
    const double factor   = 8.67373E-08;    // const. factor in mb^(-1) MeV^(-2)
    
    // Two-component Generalized Lorentzian, although almost spherical nucleus
    double GLo1 = 0.;
    double Gamma_k1 = 0.;
    double Gamma_k01 = 0.;
    double denominator1 = 0.;
    
    double GLo2 = 0.;
    double Gamma_k2 = 0.;
    double Gamma_k02 = 0.;
    double denominator2 = 0.;
    
    double GLo = 0.;
    double E1strength = 0.;
    
    // Seven parameters: par[0] = E_r1, par[1] = Gamma_r1, par[2] = sigma_r1, 
    //                   par[3] = E_r2, par[4] = Gamma_r2, par[5] = sigma_r2
    //                   par[6] = temp (common for both GLO's)
    
    // First part
    Gamma_k1 = par[1]*(pow(x[0],2.0) + (2.0*Pi*par[6])*(2.0*Pi*par[6]))/pow(par[0],2.0);
    Gamma_k01 = par[1]*((2.*Pi*par[6])*(2.*Pi*par[6]))/pow(par[0],2.0);	// Gamma_k for E_gamma=0
    denominator1 = ((pow(x[0],2.) - pow(par[0],2.))*(pow(x[0],2.) - pow(par[0],2.))) + pow(x[0],2.)*pow(Gamma_k1,2.);
    GLo1 = factor*par[2]*par[1]*((x[0]*Gamma_k1)/denominator1 + 0.7*Gamma_k01/pow(par[0],3.)); 
    
    // Second part
    Gamma_k2 = par[4]*(pow(x[0],2.0) + (2.0*Pi*par[6])*(2.0*Pi*par[6]))/pow(par[3],2.0);
    Gamma_k02 = par[4]*((2.*Pi*par[6])*(2.*Pi*par[6]))/pow(par[3],2.0);	// Gamma_k for E_gamma=0
    denominator2 = ((pow(x[0],2.) - pow(par[3],2.))*(pow(x[0],2.) - pow(par[3],2.))) + pow(x[0],2.)*pow(Gamma_k2,2.);
    GLo2 = factor*par[5]*par[4]*((x[0]*Gamma_k2)/denominator2 + 0.7*Gamma_k02/pow(par[3],3.)); 
    
    GLo = GLo1 + GLo2;
    
    E1strength = GLo;
    
    return E1strength;
}

double FitFunctionM1(double *x, double *par){
    
    const double Pi = 3.14159;
    const double factor   = 8.6737E-08;	// const. factor in mb^(-1) MeV^(-2)
    
    // Pygmy1, standard Lorentzian, M1 spin flip
    double SLo_pyg1 = 0.;
    // 3 parameters: par[0] = E_pyg1, par[1] = Gamma_pyg1, par[2] = sigma_pyg1
    double denominator_pyg1 = (pow(x[0],2.0) - pow(par[0],2.0))*(pow(x[0],2.0) - pow(par[0],2.0)) + pow(par[1],2.0)*pow(x[0],2.0);
    SLo_pyg1 = factor*par[2]*pow(par[1],2.0)*x[0]/denominator_pyg1;
    
    // Upbend function, as an exponential
    // Two parameters: const_upbend1 = par[3], const_upbend2 = par[4]
    double upbend_M1 = 0.;  
    upbend_M1 = par[3]*exp(-par[4]*x[0]);
    
    double strength_function_M1 = SLo_pyg1 + upbend_M1;
    
    return strength_function_M1;
    
}    


// Function for both E1 and M1 strength 
double FitFunctionAll(double *x, double *par){
    const double Pi = 3.14159;
    const double factor   = 8.6737E-08;	// const. factor in mb^(-1) MeV^(-2)
    
    // Two-component Generalized Lorentzian, although almost spherical nucleus
    double GLo1 = 0.;
    double Gamma_k1 = 0.;
    double Gamma_k01 = 0.;
    double denominator1 = 0.;
    
    double GLo2 = 0.;
    double Gamma_k2 = 0.;
    double Gamma_k02 = 0.;
    double denominator2 = 0.;
    
    double GLo = 0.;
    
    
    // Seven parameters: par[0] = E_r1, par[1] = Gamma_r1, par[2] = sigma_r1, 
    //                   par[3] = E_r2, par[4] = Gamma_r2, par[5] = sigma_r2
    //                   par[6] = temp (common for both GLO's)
    
    // First part, GLO
    Gamma_k1 = par[1]*(pow(x[0],2.0) + (2.0*Pi*par[6])*(2.0*Pi*par[6]))/pow(par[0],2.0);
    Gamma_k01 = par[1]*((2.*Pi*par[6])*(2.*Pi*par[6]))/pow(par[0],2.0);	// Gamma_k for E_gamma=0
    denominator1 = ((pow(x[0],2.) - pow(par[0],2.))*(pow(x[0],2.) - pow(par[0],2.))) + pow(x[0],2.)*pow(Gamma_k1,2.);
    GLo1 = factor*par[2]*par[1]*((x[0]*Gamma_k1)/denominator1 + 0.7*Gamma_k01/pow(par[0],3.)); 
    
    // Second part, GLO
    Gamma_k2 = par[4]*(pow(x[0],2.0) + (2.0*Pi*par[6])*(2.0*Pi*par[6]))/pow(par[3],2.0);
    Gamma_k02 = par[4]*((2.*Pi*par[6])*(2.*Pi*par[6]))/pow(par[3],2.0);	// Gamma_k for E_gamma=0
    denominator2 = ((pow(x[0],2.) - pow(par[3],2.))*(pow(x[0],2.) - pow(par[3],2.))) + pow(x[0],2.)*pow(Gamma_k2,2.);
    GLo2 = factor*par[5]*par[4]*((x[0]*Gamma_k2)/denominator2 + 0.7*Gamma_k02/pow(par[3],3.)); 
    
    GLo = GLo1 + GLo2;
    
    // Pygmy1, standard Lorentzian, M1 spin flip
    // 3 parameters: par[7] = E_pyg1, par[8] = Gamma_pyg1, par[9] = sigma_pyg1
    double SLo_pyg1 = 0.;
    double denominator_pyg1 = (pow(x[0],2.0) - pow(par[7],2.0))*(pow(x[0],2.0) - pow(par[7],2.0)) + pow(par[8],2.0)*pow(x[0],2.0);
    SLo_pyg1 = factor*par[9]*pow(par[8],2.0)*x[0]/denominator_pyg1;
    
    // Upbend function, as an exponential
    // Two parameters: const_upbend1 = par[10], const_upbend2 = par[11]
    double upbend_M1 = 0.;
    upbend_M1 = par[10]*exp(-par[11]*x[0]);
    
    double strength_function = GLo + SLo_pyg1 + upbend_M1;
    
    return strength_function;
   
}


void plotstrength_theory_and_all_Y_data(){
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadBorderMode(0);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(60000);
        
    // Input files, GAMMA STRENGTH FUNCTION
    ifstream strengthfile("strength.nrm");
    ifstream strengthfile_D0up("strength.nrm");
    ifstream strengthfile_D0down("strength.nrm");
    ifstream gdrfile("UnfoldedCrossSections_evaluated_at_Emax_RSFZn68.dat");
    ifstream strengthfile_low("strength.nrm");
    ifstream strengthfile_high("strength.nrm");
    ifstream transfile("transext.nrm");
    //ifstream strengthfile89Y("y89_HFB08/strength.nrm");
    //ifstream transfile89Y("y89_HFB08/transext.nrm");
    //ifstream strengthfile89Y_D0up("y89_HFB08_D0up/strength.nrm");
    //ifstream strengthfile89Y_D0down("y89_HFB08_D0down/strength.nrm");
    //ifstream strengthfile89Y_Ggup("y89_HFB08_Ggup/strength.nrm");
    //ifstream strengthfile89Y_Ggdown("y89_HFB08_Ggdown/strength.nrm");

    //ifstream gdrfile("Y89_2018_12_14_onlygn.txt");
    // From Gry, 2 Oct 2018
    //ifstream gdrfile_NEW("../gamma_n_data_89Y/unfolded_89Y.txt");
    // From Gry, 25 Jan 2019
    //ifstream gdrfile_NEW("Y89_2018_12_14_onlygn.txt");
    //ifstream Zr91gdrfile("../gamma_n_data_89Y/91Zr_Utsunomiya2008.txt");
    //ifstream Zr90Bermangdrfile("../gamma_n_data_89Y/90Zr_gn_gnp_Berman1967.txt");
    //ifstream Zr90Lepretregdrfile("../gamma_n_data_89Y/90Zr_gn_gnp_Lepretre1971.txt");
    //ifstream strength_gn_file_Lepretre("../gamma_n_data_89Y/RSF_89Y_gn_gnp_Lepretre1971.txt");
    //ifstream strength_gn_file_Berman("../gamma_n_data_89Y/RSF_89Y_gn_gnp_Berman1967.txt");
    //ifstream strength_gn_file_Young("../gamma_n_data_89Y/RSF_89Y_gn_Young1972.txt");
    //ifstream strength_gg_file_Benouaret("../gamma_n_data_89Y/RSF_89Y_gg_Benouaret2009.txt");
    
    // CALCULATIONS
    //ifstream newE1QRPAfile("../gamma_n_data_89Y/z039n051.e1");
    //ifstream newM1QRPAfile("../gamma_n_data_89Y/z039n051.m1");
    // From mail to Gry from Stephane, 15 March 2018
    //ifstream newE1QRPAfile_Gor("../gamma_n_data_89Y/z039n051_e1");
    //ifstream newM1QRPAfile_Gor("../gamma_n_data_89Y/z039n051_m1");
    // Shell-model calculations from Ronald, November 29, 2017
    // With Ex gates 3.67-7.84 MeV
    //ifstream M1_SMcalc_Exgate_file("../SM_calc_Schwengner/90Y-PN-fM1-Jav.dat");


    double strength[71],strengtherr[71],energy[100],trans[100],transerr[100];
    double strength_low[71],strengtherr_low[71];
    double strength_high[71],strengtherr_high[71];
    double strength_D0up[71],strengtherr_D0up[71];
    double strength_D0down[71],strengtherr_D0down[71];
    double error_up[71] = {0.}, error_down[71] = {0.};
    double strength_allup[71]={0.},strength_alldown[71]={0.};
    double strength_89Y[105],strengtherr_89Y[105],energy_89Y[200],transext[200];
    double strength_89Y_D0up[105],strengtherr_89Y_D0up[105],strength_89Y_D0down[105],strengtherr_89Y_D0down[105];
    double strength_89Y_Ggup[105],strengtherr_89Y_Ggup[105],strength_89Y_Ggdown[105],strengtherr_89Y_Ggdown[105];
    double error_89Y_up[105]={0.}, error_89Y_down[105]={0.};

    double e_gdr[850], gdrdata[850],gdrdataerr[850],e_gdrfit[850];
    double e_gdr_new[20], gdrdata_new[20],gdrdataerr_new[20];
    double gdrdatafit[20],gdrdataerrfit_low[20],gdrdataerrfit_high[20],gdrdataerrfit[20];
    double gdrdata_low[20], gdrdata_high[20];
    double energyerr[105] = {0.};
    double eg_Lepretre[120], gn_data_Lepretre[120], gn_dataerr_Lepretre[120];
    double eg_Berman[68], gn_data_Berman[68], gn_dataerr_Berman[68];
    double eg_Lepretre90Zr[120], gn_data_Lepretre90Zr[120], gn_dataerr_Lepretre90Zr[120];
    double eg_Berman90Zr[120], gn_data_Berman90Zr[120], gn_dataerr_Berman90Zr[120];
    double eg_Young[48], gn_data_Young[48];
    double eg_Benouaret[67], gg_data_Benouaret[67], gg_dataerr_Benouaret[67];
    double eg_91Zr[25], gdr_91Zr[25], gdrerr_91Zr[25];

    // QRPA
    double e_newqrpa[310], f_newE1qrpa[310],f_newM1qrpa[310],f_newtotqrpa[310];
    double e_newqrpa_Gor[310], f_newE1qrpa_Gor[310],f_newM1qrpa_Gor[310],f_newtotqrpa_Gor[310];
    // Shell-model
    double e_M1_strength_Exgate[105], M1_strength_Exgate[105];

    int i,j;
    double a0 =  -0.7715;
    double a1 =   0.1267;
    double a0_89Y =  -0.9411;
    double a1_89Y =   0.1228;
    double a0_68Zn = -0.8360;
    double a1_68Zn = 0.1280;

    double eg_err[900] = {0.};
    double pi = 3.14159265359;

    
    string line;
    line.resize(256);
    
    double x,y,z,d1,d2;	
    const double factor = 8.674E-08;	// const. factor in mb^(-1) MeV^(-2)


    // Read Oslo data, 90Y
    // Recommended norm
    i=0;
    while(strengthfile){
        strengthfile >> x;
        if(i<70){
            strength[i] = x;
            energy[i] = a0 + (a1*i);
            cout << i << " " << energy[i] << " " << strength[i] << endl;
        }	
        else{
            strengtherr[i-70] = x;
        }
        i++;
    }
    strengthfile.close();

    // Read transmission coefficient and transform to gSF to check extrapolations
    for(i=0;i<70;i++){
        transfile >> x;
        energy[i] = a0 + (a1*i);
        trans[i] = x/(2.*pi*pow(energy[i],3.));
        transerr[i] = 0.1*trans[i];

    }
    transfile.close();
    
    // Ggdown
    i=0;
    while(strengthfile_low){
        strengthfile_low >> x;
        if(i<70){
            strength_low[i] = x;
        }	
        else{strengtherr_low[i-70] = x;}
        i++;
    }
    strengthfile_low.close();
    
    // Ggup
    i=0;
    while(strengthfile_high){
        strengthfile_high >> x;
        if(i<70){
            strength_high[i] = x;
        }	
        else{strengtherr_high[i-70] = x;}
        i++;
    }
    strengthfile_high.close();

    // D0up
    i=0;
    while(strengthfile_D0up){
        strengthfile_D0up >> x;
        if(i<70){
            strength_D0up[i] = x;
        }	
        else{strengtherr_D0up[i-70] = x;}
        i++;
    }
    strengthfile_D0up.close();

    // D0down
    i=0;
    while(strengthfile_D0down){
        strengthfile_D0down >> x;
        if(i<70){
            strength_D0down[i] = x;
        }   
        else{strengtherr_D0down[i-70] = x;}
        i++;
    }
    strengthfile_D0down.close();
    
    // Estimate errors, approx. 1sigma, 90Y
    for(i=0;i<70;i++){
        // Standard deviation from recommended
        if(strength[i]>0.) {
            error_up[i] = strength[i] * sqrt( pow(((strength_D0up[i]-strength[i])/strength[i]),2.) 
                                             + pow(((strength_high[i]-strength[i])/strength[i]),2.)
                                             + pow((strengtherr[i]/strength[i]),2.)
                                             );
            error_down[i] = strength[i] * sqrt( pow(((strength_D0down[i]-strength[i])/strength[i]),2.) 
                                             + pow(((strength_low[i]-strength[i])/strength[i]),2.) 
                                             + pow((strengtherr[i]/strength[i]),2.)
                                             );
        }        
        else {
            error_up[i] = strengtherr[i];
            error_down[i] = strengtherr[i];
        }  
        strength_allup[i] = strength[i] + error_up[i];
        strength_alldown[i] = strength[i] - error_down[i];
    }

    

    // Read Oslo data, 89Y  
    // Recommended normalization


    //i=0;
    //while(strengthfile89Y){
    //    strengthfile89Y >> x;
    //    if(i<99){
    //        strength_89Y[i] = x;
    //        energy_89Y[i] = a0_89Y + (a1_89Y*i);
    //    }
    //    else{strengtherr_89Y[i-99] = x;}
    //    i++;
    //}
    //strengthfile89Y.close();
//
    //// Read transmission coefficient and transform to gSF to check extrapolations
    //for(i=0;i<162;i++){
    //    transfile89Y >> x;
    //    energy_89Y[i] = a0_89Y + (a1_89Y*i);
    //    transext[i] = x/(2.*pi*pow(energy_89Y[i],3.));
    //}
    //transfile89Y.close();
    //
    //// D0up
    //i=0;
    //while(strengthfile89Y_D0up){
    //    strengthfile89Y_D0up >> x;
    //    if(i<99){
    //        strength_89Y_D0up[i] = x;
    //    }
    //    else{strengtherr_89Y_D0up[i-99] = x;}
    //    i++;
    //}
    //strengthfile89Y_D0up.close();
    //
    //// D0down
    //i=0;
    //while(strengthfile89Y_D0down){
    //    strengthfile89Y_D0down >> x;
    //    if(i<99){
    //        strength_89Y_D0down[i] = x;
    //    }
    //    else{strengtherr_89Y_D0down[i-99] = x;}
    //    i++;
    //}
    //strengthfile89Y_D0down.close();
//
    //// Ggup
    //i=0;
    //while(strengthfile89Y_Ggup){
    //    strengthfile89Y_Ggup >> x;
    //    if(i<99){
    //        strength_89Y_Ggup[i] = x;
    //    }
    //    else{strengtherr_89Y_Ggup[i-99] = x;}
    //    i++;
    //}
    //strengthfile89Y_Ggup.close();
//
    //// Ggdown
    //i=0;
    //while(strengthfile89Y_Ggdown){
    //    strengthfile89Y_Ggdown >> x;
    //    if(i<99){
    //        strength_89Y_Ggdown[i] = x;
    //    }
    //    else{strengtherr_89Y_Ggdown[i-99] = x;}
    //    i++;
    //}
    //strengthfile89Y_Ggdown.close();
    
    // Estimate errors, approx. 1sigma, 89Y
    //for(i=0;i<99;i++){
    //    // Standard deviation from recommended
    //    if(strength_89Y[i]>0.) {
    //        error_89Y_up[i] = strength_89Y[i] * sqrt( pow(((strength_89Y_Ggup[i]-strength_89Y[i])/strength_89Y[i]),2.) 
    //                                         + pow(((strength_89Y_D0up[i]-strength_89Y[i])/strength_89Y[i]),2.)
    //                                         + pow((strengtherr_89Y[i]/strength_89Y[i]),2.)
    //                                         );
    //        error_89Y_down[i] = strength_89Y[i] * sqrt( pow(((strength_89Y_D0down[i]-strength_89Y[i])/strength_89Y[i]),2.) 
    //                                           + pow(((strength_89Y_Ggdown[i]-strength_89Y[i])/strength_89Y[i]),2.) 
    //                                           + pow((strengtherr_89Y[i]/strength_89Y[i]),2.)
    //                                           );
    //    }
    //    else {
    //        error_89Y_up[i] = strengtherr_89Y[i];
    //        error_89Y_down[i] = strengtherr_89Y[i];
    //    }  
    //    //cout << i << " " << strength_89Y[i] << " " << strengtherr_89Y[i] << " " << error_89Y_down[i] << " " << error_89Y_up[i] << endl;
    //}
    
    

    

    // Read Japan data
    i=0;
    while(gdrfile){
        gdrfile >> x >> y;
        e_gdr[i] = x; // in MeV
        gdrdata[i] = factor*y/x;
        gdrdataerr[i] = 0.045*factor*y/x; //  4.5% error corresponds to 1sigma according to Gry (email 19 Feb 2018)
        //cout << i << " " << x << endl;
        i++;
    }
    gdrfile.close();

    // NEW Japan data, file from Gry on Oct 2, 2018
    /*getline(gdrfile_NEW,line); // skip first line with text

    i=0;
    while(gdrfile_NEW){
        gdrfile_NEW >> x >> y  >> z;
        e_gdr_new[i] = x; // in MeV
        gdrdata_new[i] = factor*y/x;
        gdrdataerr_new[i] = factor*z/x; 
        gdrdata_low[i] =  (factor*y/x) - (factor*z/x);
        gdrdata_high[i] =  (factor*y/x) + (factor*z/x);
        i++;
    }
    gdrfile_NEW.close();
    */

    // NEW Japan data, file from Gry on Jan 25, 2019
    i=0;
    while(gdrfile_NEW){
        gdrfile_NEW >> x >> y  >> z;
        e_gdr_new[i] = x; // in MeV
        gdrdata_new[i] = factor*y/x;
        gdrdataerr_new[i] = factor*z/x; 
        gdrdata_low[i] =  (factor*y/x) - (factor*z/x);
        gdrdata_high[i] =  (factor*y/x) + (factor*z/x);
        //cout << i << " " << e_gdr_new[i] << " " << gdrdata_new[i] << endl;
        if(i>3){
            e_gdrfit[i-4] = e_gdr_new[i];
            gdrdatafit[i-4] = gdrdata_new[i];
            gdrdataerrfit[i-4] = gdrdataerr_new[i];
            gdrdataerrfit_low[i-4] = gdrdataerr_new[i];
            gdrdataerrfit_high[i-4] =gdrdataerr_new[i];
            cout << i-4 << " " << e_gdrfit[i-4] << " " << gdrdatafit[i-4] << endl;
        } 
        i++;
    }
    gdrfile_NEW.close();

    // 91Zr Japan data
    for(i=0;i<13;i++){    // Skip text
        getline(Zr91gdrfile,line);
    }
    
    i=0;
    while(Zr91gdrfile){
        Zr91gdrfile >> x >> y >> z;
        eg_91Zr[i] = x; // in MeV
        gdr_91Zr[i] = factor*y/x;
        gdrerr_91Zr[i] = factor*z/x; // 
        //cout << eg_91Zr[i] << " " << gdr_91Zr[i] << " " << gdrerr_91Zr[i] << endl;
        i++;
    }
    Zr91gdrfile.close();

    
    // More (GAMMA,N) DATA
	//getline(strength_gn_file_Lepretre,line);
	//getline(strength_gn_file_Lepretre,line);
	//getline(strength_gn_file_Lepretre,line);
	//getline(strength_gn_file_Lepretre,line);
	//getline(strength_gn_file_Lepretre,line);
	//getline(strength_gn_file_Lepretre,line);
	//getline(strength_gn_file_Lepretre,line);
	//getline(strength_gn_file_Lepretre,line);
	//getline(strength_gn_file_Lepretre,line);
	//i=0;
    //while(strength_gn_file_Lepretre){
    //    strength_gn_file_Lepretre >> x >> y >> z;
    //    eg_Lepretre[i]         = x;
    //    gn_data_Lepretre[i]    = y;
    //    gn_dataerr_Lepretre[i] = z;
    //    i++;
    //}
    //strength_gn_file_Lepretre.close();
//
//
	//getline(strength_gn_file_Berman,line);
	//getline(strength_gn_file_Berman,line);
	//getline(strength_gn_file_Berman,line);
    //i=0;
    //while(strength_gn_file_Berman){
    //    strength_gn_file_Berman >> x >> y >> z;
    //    eg_Berman[i]         = x;
    //    gn_data_Berman[i]    = y;
    //    gn_dataerr_Berman[i] = z;
    //    i++;
    //}
    //strength_gn_file_Berman.close();
    //
    //i=0;
    //while(strength_gn_file_Young){
    //    strength_gn_file_Young >> x >> y;
    //    eg_Young[i]         = x;
    //    gn_data_Young[i]    = y;
    //    i++;
    //}
    //strength_gn_file_Young.close();
    //
    //i=0;
    //while(strength_gg_file_Benouaret){
    //    strength_gg_file_Benouaret >> x >> y >> z;
    //    eg_Benouaret[i]         = x;
    //    gg_data_Benouaret[i]    = y;
    //    gg_dataerr_Benouaret[i] = z;
    //    i++;
    //}
    //strength_gg_file_Benouaret.close();
//
    //// 90Zr Lepretre data
    //for(i=0;i<10;i++){    // Skip text
    //    getline(Zr90Lepretregdrfile,line);
    //}
    //
    //i=0;
    //while(Zr90Lepretregdrfile){
    //    Zr90Lepretregdrfile >> x >> y >> z;
    //    eg_Lepretre90Zr[i] = x; // in MeV
    //    gn_data_Lepretre90Zr[i] = factor*y/x;
    //    gn_dataerr_Lepretre90Zr[i] = factor*z/x; // 
    //    //cout << eg_91Zr[i] << " " << gdr_91Zr[i] << " " << gdrerr_91Zr[i] << endl;
    //    i++;
    //}
    //Zr90Lepretregdrfile.close();
//
    //// 90Zr Berman data
   f//or(i=0;i<12;i++){    // Skip text
    //    getline(Zr90Bermangdrfile,line);
    //}
    //
    //i=0;
    //while(Zr90Bermangdrfile){
    //    Zr90Bermangdrfile >> x >> y >> z;
    //    eg_Berman90Zr[i] = x; // in MeV
    //    gn_data_Berman90Zr[i] = factor*y/x;
    //    gn_dataerr_Berman90Zr[i] = factor*z/x; // 
    //    //cout << eg_91Zr[i] << " " << gdr_91Zr[i] << " " << gdrerr_91Zr[i] << endl;
    //    i++;
    //}
    //Zr90Bermangdrfile.close();
    
    // E1 Gamma strengths from RIPL-2 (Kopecky evalutation)
    double e_74Ge[1] = {7.1}; // MeV
    double f_74Ge[1] = {3.44E-08}; // MeV^-3
    double ferr_74Ge[1] = {1.15E-08}; // MeV^-3

    double e_91Zr[1] = {7.2};
    double f_91Zr[1] = {3.14E-08};
    double ferr_91Zr[1] = {1.18E-08};

    double e_94Nb[1] = {6.5};
    double f_94Nb[1] = {2.84E-08};
    double ferr_94Nb[] = {0.7E-08};
    
    // New caclulations from Stephane Goriely (email 10 May 2017)
    // E1
    getline(newE1QRPAfile,line);
    getline(newE1QRPAfile,line);
    i=0;
    while(newE1QRPAfile){
        newE1QRPAfile >> x >> y;
        //cout << i << " " << x << " " << y << endl;
        e_newqrpa[i] = x;
        f_newE1qrpa[i] = factor*y;
        i++;
    }
    newE1QRPAfile.close();
    // M1
    getline(newM1QRPAfile,line);
    getline(newM1QRPAfile,line);
    i=0;
    while(newM1QRPAfile){
        newM1QRPAfile >> x >> y;
        //cout << i << " " << x << " " << y << endl;
        f_newM1qrpa[i] = factor*y;
        f_newtotqrpa[i] = f_newE1qrpa[i] + factor*y;
        i++;
    }
    newM1QRPAfile.close();
    
    i=0;
    while(M1_SMcalc_Exgate_file){
        M1_SMcalc_Exgate_file >> x >> y;
        e_M1_strength_Exgate[i] = x;
        M1_strength_Exgate[i] = y;
        i++;
    }
    M1_SMcalc_Exgate_file.close();
    
    // New caclulations from Stephane Goriely (email to Gry 15 March 2018)
    // E1
    getline(newE1QRPAfile_Gor,line);
    getline(newE1QRPAfile_Gor,line);
    i=0;
    while(newE1QRPAfile_Gor){
        newE1QRPAfile_Gor >> x >> y;
        //cout << i << " " << x << " " << y << endl;
        e_newqrpa_Gor[i] = x;
        f_newE1qrpa_Gor[i] = factor*y;
        i++;
    }
    newE1QRPAfile_Gor.close();
    // M1
    getline(newM1QRPAfile_Gor,line);
    getline(newM1QRPAfile_Gor,line);
    i=0;
    while(newM1QRPAfile_Gor){
        newM1QRPAfile_Gor >> x >> y;
        //cout << i << " " << x << " " << y << endl;
        f_newM1qrpa_Gor[i] = factor*y;
        f_newtotqrpa_Gor[i] = f_newE1qrpa_Gor[i] + factor*y;
        i++;
    }
    newM1QRPAfile_Gor.close();
    
    // Shell-model calculations
    i=0;
    while(M1_SMcalc_Exgate_file){
        M1_SMcalc_Exgate_file >> x >> y;
        e_M1_strength_Exgate[i] = x;
        M1_strength_Exgate[i] = y;
        i++;
    }
    M1_SMcalc_Exgate_file.close();
    
    

    
    
    // PLOT THINGS
    // Make graphs
    TGraphAsymmErrors *strengthexp = new TGraphAsymmErrors(69,energy,strength,energyerr,energyerr,strengtherr,strengtherr);
    TGraph *strengthexp_D0up = new TGraph(69,energy,strength_D0up);
    TGraph *strengthexp_D0down = new TGraph(69,energy,strength_D0down);
    TGraph *strengthexp_Ggup = new TGraph(69,energy,strength_high);
    TGraph *strengthexp_Ggdown = new TGraph(69,energy,strength_low);
    TGraph *strengthexp_allup = new TGraph(69,energy,strength_allup);
    TGraph *strengthexp_alldown = new TGraph(69,energy,strength_alldown);

    TGraphErrors *transexp = new TGraphErrors(70,energy,trans,energyerr,transerr); // The whole thing
    //TGraphErrors *transexp = new TGraphErrors(21,energy,trans,energyerr,transerr); // Only low-energy extrapolation
    TGraphErrors *strengthexp_rec = new TGraphErrors(69,energy,strength,energyerr,strengtherr);
    TGraphErrors *strengthexp_low = new TGraphErrors(69,energy,strength_low,energyerr,strengtherr_low);
    TGraphErrors *strengthexp_high = new TGraphErrors(69,energy,strength_high,energyerr,strengtherr_high);
    TGraphAsymmErrors *strengthexp_89Y = new TGraphAsymmErrors(98,energy_89Y,strength_89Y,energyerr,energyerr,error_89Y_down,error_89Y_up);

    TGraphErrors *gdrstrength = new TGraphErrors(840,e_gdr,gdrdata,eg_err,gdrdataerr);
    TGraphErrors *gdrstrengthfit = new TGraphErrors(15,e_gdrfit,gdrdatafit,eg_err,gdrdataerrfit);

    TGraphErrors *gdrstrength_new = new TGraphErrors(17,e_gdr_new,gdrdata_new,eg_err,gdrdataerr_new);
    TGraph *gdrstrength_low = new TGraph(17,e_gdr_new,gdrdata_low);
    TGraph *gdrstrength_high = new TGraph(17,e_gdr_new,gdrdata_high);
 
    TGraphErrors *gammaexp_89Y_Lepretre = new TGraphErrors(109,eg_Lepretre,gn_data_Lepretre,eg_err,gn_dataerr_Lepretre);
    TGraphErrors *gammaexp_89Y_Berman = new TGraphErrors(64,eg_Berman,gn_data_Berman,eg_err,gn_dataerr_Berman);
    TGraph *gammaexp_89Y_Young = new TGraph(47,eg_Young,gn_data_Young);
    TGraphErrors *gammaexp_89Y_Benouaret = new TGraphErrors(66,eg_Benouaret,gg_data_Benouaret,eg_err,gg_dataerr_Benouaret);

    
    TGraphErrors *gdrstrength91Zr = new TGraphErrors(24,eg_91Zr,gdr_91Zr,eg_err,gdrerr_91Zr);
    TGraphErrors *gammaexp_90Zr_Lepretre = new TGraphErrors(119,eg_Lepretre90Zr,gn_data_Lepretre90Zr,eg_err,gn_dataerr_Lepretre90Zr);
    TGraphErrors *gammaexp_90Zr_Berman = new TGraphErrors(90,eg_Berman90Zr,gn_data_Berman90Zr,eg_err,gn_dataerr_Berman90Zr);
    
    TGraphErrors *strengthexp_74Ge = new TGraphErrors(1,e_74Ge,f_74Ge,energyerr,ferr_74Ge);
    TGraphErrors *strengthexp_91Zr = new TGraphErrors(1,e_91Zr,f_91Zr,energyerr,ferr_91Zr);
    TGraphErrors *strengthexp_94Nb = new TGraphErrors(1,e_94Nb,f_94Nb,energyerr,ferr_94Nb);
    
    TGraph *trans_89Y = new TGraph(106,energy_89Y,transext);
   
    TGraph *newTotqrpa_graph = new TGraph(280,e_newqrpa,f_newtotqrpa);
    TGraph *newM1qrpa_graph = new TGraph(280,e_newqrpa,f_newM1qrpa);
    TGraph *newE1qrpa_graph = new TGraph(280,e_newqrpa,f_newE1qrpa);
    TGraph *newTotqrpa_Gor_graph = new TGraph(280,e_newqrpa_Gor,f_newtotqrpa_Gor);
    TGraph *newM1qrpa_Gor_graph = new TGraph(280,e_newqrpa_Gor,f_newM1qrpa_Gor);
    TGraph *newE1qrpa_Gor_graph = new TGraph(280,e_newqrpa_Gor,f_newE1qrpa_Gor);
    TGraph *M1_SM_strength_Exgate = new TGraphErrors(37,e_M1_strength_Exgate,M1_strength_Exgate);
    
    // Make a combined TGraphErrors for OCL data and GDR data 
    TMultiGraph *OCL_and_GDR_data = new TMultiGraph();
    //OCL_and_GDR_data->Add(transexp,"P");
    //OCL_and_GDR_data->Add(strengthexp_low,"P");
    OCL_and_GDR_data->Add(strengthexp_rec,"P");
    //OCL_and_GDR_data->Add(strengthexp_high,"P");
    OCL_and_GDR_data->Add(gdrstrengthfit,"P");
    

    TCanvas *c1 = new TCanvas("c1","Gamma-ray strength function",1200,500);
    c1->Divide(2,1,0,0);
    TH2F *h = new TH2F("h"," ",100,0.0,   19.155,100,6.480e-10,2.167e-06);
    TH2F *h2 = new TH2F("h2"," ",100,0.0,   19.155,100,6.480e-10,2.167e-06);

    // First panel
    c1->cd(1);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitle("#gamma-ray energy E_{#gamma} (MeV)");
    h->GetXaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->SetTitleFont(42);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->SetTitleSize(0.043);
    h->GetYaxis()->CenterTitle();
    h->GetYaxis()->SetTitleOffset(1.3);
    h->GetYaxis()->SetTitle("#gamma-ray strength function (MeV^{-3})");
    h->GetYaxis()->SetTitleFont(42);
    h->GetYaxis()->SetLabelFont(42);
    h->GetYaxis()->SetTitleSize(0.043);
    h->Draw();
    
    transexp->SetMarkerStyle(25);
    //transexp->Draw("P");
    
    // Error due to <Gg>
    strengthexp_low->SetFillStyle(1001);
    strengthexp_low->SetFillColor(kMagenta-10);
    strengthexp_low->SetLineWidth(402);
    strengthexp_low->SetLineColor(kMagenta-10);
    strengthexp_low->DrawGraph(51,&strengthexp_low->GetX()[18],&strengthexp_low->GetY()[18],"L");
    
    strengthexp_high->SetFillStyle(1001);
    strengthexp_high->SetFillColor(kMagenta-10);
    strengthexp_high->SetLineWidth(-402);
    strengthexp_high->SetLineColor(kMagenta-10);
    strengthexp_high->DrawGraph(51,&strengthexp_high->GetX()[18],&strengthexp_high->GetY()[18],"L");

    // Error due to D0
    strengthexp_D0down->SetFillStyle(1001);
    strengthexp_D0down->SetFillColor(kMagenta-3);
    strengthexp_D0down->SetLineWidth(101);
    strengthexp_D0down->SetLineColor(kMagenta-3);
    //strengthexp_D0down->SetMarkerStyle(2);
    strengthexp_D0down->DrawGraph(51,&strengthexp_D0down->GetX()[18],&strengthexp_D0down->GetY()[18],"LP");

    strengthexp_D0up->SetFillStyle(1001);
    strengthexp_D0up->SetFillColor(kMagenta-3);
    strengthexp_D0up->SetLineWidth(-101);
    strengthexp_D0up->SetLineColor(kMagenta-3);
    strengthexp_D0up->DrawGraph(51,&strengthexp_D0up->GetX()[18],&strengthexp_D0up->GetY()[18],"L");

    // Total errors assuming independent uncertainties
    strengthexp_allup->SetLineColor(kMagenta);
    strengthexp_allup->SetLineStyle(7);
    strengthexp_allup->DrawGraph(51,&strengthexp_allup->GetX()[18],&strengthexp_allup->GetY()[18],"L");

    strengthexp_alldown->SetLineColor(kMagenta);
    strengthexp_alldown->SetLineStyle(7);
    strengthexp_alldown->DrawGraph(51,&strengthexp_alldown->GetX()[18],&strengthexp_alldown->GetY()[18],"L");

    // 89Y (not used)
    strengthexp_89Y->SetMarkerStyle(22);
    strengthexp_89Y->SetMarkerSize(0.8);
    strengthexp_89Y->SetMarkerColor(kPink+1);
    strengthexp_89Y->SetLineColor(kPink+1);
    //strengthexp_89Y->Draw("P");
 
    strengthexp->SetMarkerStyle(21);
    strengthexp->SetMarkerSize(0.6);
    strengthexp->Draw("P");
    

    trans_89Y->SetMarkerStyle(25);
    //trans_89Y->Draw("P");


    gdrstrength->SetMarkerStyle(2);
    gdrstrength->SetMarkerSize(1.);
    //gdrstrength->SetMarkerColor(kAzure+7);
    //gdrstrength->SetLineColor(kAzure+7);
    //gdrstrength->SetLineWidth(3);
    //gdrstrength->Draw("L");
    //gdrstrength->DrawGraph(700,&gdrstrength->GetX()[70],&gdrstrength->GetY()[70],"P");
   
    gdrstrengthfit->SetMarkerStyle(2);
    gdrstrengthfit->SetMarkerSize(1.5);
    gdrstrengthfit->SetMarkerColor(kBlue);
    gdrstrengthfit->SetLineColor(kBlue);
    //gdrstrengthfit->Draw("P");

    gdrstrength_new->SetMarkerStyle(2);
    gdrstrength_new->SetLineWidth(3);
    gdrstrength_new->SetLineColor(kCyan+1);
    gdrstrength_new->DrawGraph(12,&gdrstrength_new->GetX()[3],&gdrstrength_new->GetY()[3],"LP");

    gdrstrength_low->DrawGraph(12,&gdrstrength_low->GetX()[3],&gdrstrength_low->GetY()[3],"L");
    gdrstrength_high->DrawGraph(12,&gdrstrength_high->GetX()[3],&gdrstrength_high->GetY()[3],"L");

    
    gammaexp_89Y_Young->SetMarkerStyle(22);
    gammaexp_89Y_Young->SetMarkerColor(kCyan+1);
    gammaexp_89Y_Young->SetLineColor(kCyan+1);
    gammaexp_89Y_Young->SetMarkerSize(0.8);
    //gammaexp_89Y_Young->Draw("P");

    gammaexp_89Y_Benouaret->SetMarkerStyle(24);
    gammaexp_89Y_Benouaret->SetMarkerSize(0.6);
    gammaexp_89Y_Benouaret->SetMarkerColor(kViolet+8);
    gammaexp_89Y_Benouaret->SetLineColor(kViolet+8);
    //gammaexp_89Y_Benouaret->Draw("P");


    strengthexp_74Ge->SetMarkerStyle(34);
    strengthexp_74Ge->SetMarkerSize(1.5);
    strengthexp_74Ge->SetLineColor(kPink+10);
    strengthexp_74Ge->SetMarkerColor(kPink+10);
    //strengthexp_74Ge->Draw("P");
    
    strengthexp_91Zr->SetMarkerStyle(34);
    strengthexp_91Zr->SetMarkerSize(1.5);
    strengthexp_91Zr->SetLineColor(kViolet-1);
    strengthexp_91Zr->SetMarkerColor(kViolet-1);
    //strengthexp_91Zr->Draw("P");

    strengthexp_94Nb->SetMarkerStyle(34);
    strengthexp_94Nb->SetMarkerSize(1.5);
    strengthexp_94Nb->SetLineColor(kCyan+2);
    strengthexp_94Nb->SetMarkerColor(kCyan+2);
    //strengthexp_94Nb->Draw("P");

    gammaexp_89Y_Lepretre->SetMarkerStyle(2);
    gammaexp_89Y_Lepretre->SetMarkerColor(kAzure+7);
    gammaexp_89Y_Lepretre->SetLineColor(kAzure+7);
    gammaexp_89Y_Lepretre->SetMarkerSize(1.);
    //gammaexp_89Y_Lepretre->Draw("P");
    
    gammaexp_89Y_Berman->SetMarkerStyle(28);
    gammaexp_89Y_Berman->SetMarkerSize(1.1);
    gammaexp_89Y_Berman->SetMarkerColor(kMagenta+1);
    gammaexp_89Y_Berman->SetLineColor(kMagenta+1);
    //gammaexp_89Y_Berman->Draw("P");
    
    gdrstrength91Zr->SetMarkerStyle(27);
    gdrstrength91Zr->SetMarkerSize(1.1);
    //gdrstrength91Zr->Draw("P");

    gammaexp_90Zr_Lepretre->SetMarkerStyle(24);
    gammaexp_90Zr_Lepretre->SetMarkerColor(kGreen+1);
    gammaexp_90Zr_Lepretre->SetLineColor(kGreen+1);
    //gammaexp_90Zr_Lepretre->Draw("P");

    gammaexp_90Zr_Berman->SetMarkerStyle(28);
    //gammaexp_90Zr_Berman->Draw("P");
    
    newTotqrpa_graph->SetLineWidth(1);
    newTotqrpa_graph->SetLineColor(kMagenta);
    //newTotqrpa_graph->Draw("L");
    
    newM1qrpa_graph->SetLineStyle(3);
    //newM1qrpa_graph->Draw("L same");
    
    newE1qrpa_graph->SetLineStyle(5);
    //newE1qrpa_graph->Draw("L same");
    
    //newTotqrpa_Gor_graph->SetLineColor(kMagenta+1);
    newTotqrpa_Gor_graph->Draw("L");

    newM1qrpa_Gor_graph->SetLineStyle(3);
    newM1qrpa_Gor_graph->Draw("L same");
    
    newE1qrpa_Gor_graph->SetLineStyle(5);
    newE1qrpa_Gor_graph->Draw("L same");

    M1_SM_strength_Exgate->SetLineColor(kBlue);
    M1_SM_strength_Exgate->SetLineWidth(2);
    M1_SM_strength_Exgate->SetLineStyle(2);
    M1_SM_strength_Exgate->Draw("L");


    c1->Update();

    TLegend *leg2 = new TLegend(0.22,0.65,0.43,0.93);
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.037);
    leg2->SetFillColor(0);
    leg2->AddEntry(gdrstrength," ^{89}Y(#gamma,n), this work","LP");
    //leg2->AddEntry(strengthexp_89Y," (p,p'#gamma)^{89}Y","p"); // low = 101 meV, high = 302 meV
    leg2->AddEntry(strengthexp," ^{90}Y, this work ","P"); // <Gg> = 168 meV
    leg2->AddEntry(strengthexp_D0up," ^{90}Y, uncertainty D_{0}","f"); // 
    leg2->AddEntry(strengthexp_high," ^{90}Y, uncertainty #LT#Gamma_{#gamma}#GT","f"); // 
    leg2->AddEntry(strengthexp_allup," ^{90}Y, total uncertainty ","L"); // 
    //leg2->AddEntry(newTotqrpa_graph," QRPA E1 + M1 (previous) ","L");
    leg2->Draw();
 
    TLegend *leg3 = new TLegend(0.7,0.25,0.93,0.45);
    leg3->SetBorderSize(0);
    leg3->SetTextFont(42);
    leg3->SetTextSize(0.037);
    leg3->SetFillColor(0);
     leg3->AddEntry(M1_SM_strength_Exgate," M1 shell model ","L");
    leg3->AddEntry(newE1qrpa_Gor_graph," QRPA E1 ","L");
    leg3->AddEntry(newM1qrpa_Gor_graph," QRPA M1 ","L");
    leg3->AddEntry(newTotqrpa_Gor_graph," QRPA E1 + M1 ","L");
    leg3->Draw();
   
    TLatex t;
    t.SetTextSize(0.04);
    t.SetTextFont(42);
    t.DrawLatex(1.,1.2E-06,"(a)");

    
    // Next panel
    c1->cd(2);
    gPad->SetLogy();
    gPad->SetLeftMargin(0.07);
    gPad->SetRightMargin(0.03);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.05);

    h2->GetXaxis()->CenterTitle();
    h2->GetXaxis()->SetTitle("#gamma-ray energy E_{#gamma} (MeV)");
    h2->GetXaxis()->SetTitleOffset(1.1);
    h2->GetXaxis()->SetTitleFont(42);
    h2->GetXaxis()->SetLabelFont(42);
    h2->GetXaxis()->SetTitleSize(0.045);
    h2->GetYaxis()->SetLabelFont(42);
    h2->Draw();

    strengthexp_low->DrawGraph(51,&strengthexp_low->GetX()[18],&strengthexp_low->GetY()[18],"L");
    strengthexp_high->DrawGraph(51,&strengthexp_high->GetX()[18],&strengthexp_high->GetY()[18],"L");
    strengthexp_D0down->DrawGraph(51,&strengthexp_D0down->GetX()[18],&strengthexp_D0down->GetY()[18],"LP");
    strengthexp_D0up->DrawGraph(51,&strengthexp_D0up->GetX()[18],&strengthexp_D0up->GetY()[18],"L");
    strengthexp_allup->DrawGraph(51,&strengthexp_allup->GetX()[18],&strengthexp_allup->GetY()[18],"L");
    strengthexp_alldown->DrawGraph(51,&strengthexp_alldown->GetX()[18],&strengthexp_alldown->GetY()[18],"L");

    //strengthexp_89Y->Draw("P");
    strengthexp->Draw("P");

    //transexp->Draw("P");

    strengthexp_rec->SetMarkerStyle(2);
    strengthexp_rec->SetMarkerColor(kPink+10);
    strengthexp_rec->SetLineColor(kPink+10);
    //strengthexp_rec->Draw("P");
    
    //gdrstrength->DrawGraph(700,&gdrstrength->GetX()[70],&gdrstrength->GetY()[70],"P");
    
    // NEW DATA FROM GRY
    gdrstrength_new->DrawGraph(12,&gdrstrength_new->GetX()[3],&gdrstrength_new->GetY()[3],"LP");
    gdrstrength_low->DrawGraph(12,&gdrstrength_low->GetX()[3],&gdrstrength_low->GetY()[3],"L");
    gdrstrength_high->DrawGraph(12,&gdrstrength_high->GetX()[3],&gdrstrength_high->GetY()[3],"L");
    
    //gdrstrengthfit->Draw("P");

    // START THE FITS
    // General constants
    const double Pi = 3.14159;
    
    // START VALUES for the GDR parameters. 
    double E_r1 = 1.71175e+01; // Start value in MeV
    double Gamma_r1 = 1.83983e+00; // Start value in MeV
    double sigma_r1 = 1.42160e+02; // Start value in mb
    
    double E_r2 = 1.35878e+01; // Start value in MeV
    double Gamma_r2 = 4.01777e+00; // Start value in MeV
    double sigma_r2 = 1.24612e+02; // Start value in mb
    
    // START VALUES pygmy resonance 1 
    double E_pyg1 = 6.5; // Pygmy centroid (MeV)
    double Gamma_pyg1 = 0.5; // Width (MeV)
    double sigma_pyg1 = 0.68; // peak cross section (mb)
    
    // START VALUE temperature
    double temp = 7.51451e-01; // [MeV] 
    
    // Start value for upbend constants
    double constup1 = 1.E-08;
    double constup2 = 1.5;
    
    TF1 *fit_strength_E1 = new TF1("fit_strength_E1",FitFunctionE1,1.5,18.6,7); // GDR (7 parameters)
    //Vector for 7 start parameters (GDR)
    double par_array[7] = {16.9,2.,236.,14.5586,3.53,50.,0.828014};
    //                      p0  p1  p2    p3     p4   p5   p6  
    fit_strength_E1->SetParameters(par_array);
    // Possible to set Limits on the temperature parameter to ease the fit
    //fit_strength_E1->FixParameter(6,0.95); // From CT interpolation, 90Y
 
    // Some constraints are needed for the upper-limit fit only
    //fit_strength_E1->SetParLimits(0,16.,18.);
    //fit_strength_E1->SetParLimits(1,1.5,5.);
    //fit_strength_E1->SetParLimits(3,13.,15.);
    //fit_strength_E1->SetParLimits(4,1.5,5.);
    
   
    cout << " Fit of only E1 part:" << endl;    
    //gdrstrengthfit->Fit(fit_strength_E1,"RM+");
    OCL_and_GDR_data->Fit(fit_strength_E1,"RM+");
    
    fit_strength_E1->SetLineColor(kBlack);
    fit_strength_E1->SetLineWidth(1);
    fit_strength_E1->SetLineStyle(1);
    //fit_strength_E1->Draw("L same");
    
    // Now we have the parameters for the GDR. We can fix them for the next fit, 
    // where the M1 contribution is determined
    E_r1 = fit_strength_E1->GetParameter(0);
    Gamma_r1 = fit_strength_E1->GetParameter(1);
    sigma_r1 = fit_strength_E1->GetParameter(2);
    E_r2 = fit_strength_E1->GetParameter(3);
    Gamma_r2 = fit_strength_E1->GetParameter(4);
    sigma_r2 = fit_strength_E1->GetParameter(5);
    temp = fit_strength_E1->GetParameter(6);    


    // FIT OF GDR + PYGMY + UPBEND
    //Vector for 12 start parameters (GDR, pygmy, and upbend)
    double par2_array[12] = {E_r1,Gamma_r1,sigma_r1,E_r2,Gamma_r2,sigma_r2,temp,E_pyg1,Gamma_pyg1,sigma_pyg1,constup1,constup2};
    //                        p0    p1      p2       p3      p4      p5     p6    p7      p8        p9        p10       p11      
    
    
    TF1 *fit_strength_all = new TF1("fit_strength_all",FitFunctionAll,1.5,18.6,12); // GDR, pygmy and upbend (12 parameters, but only 11 free parameters)
    fit_strength_all->SetParameters(par2_array);
    
    // Fix first 7 parameters from the previous fit
    /*fit_strength_all->FixParameter(0,E_r1); 
    fit_strength_all->FixParameter(1,Gamma_r1); 
    fit_strength_all->FixParameter(2,sigma_r1); 
    fit_strength_all->FixParameter(3,E_r2); 
    fit_strength_all->FixParameter(4,Gamma_r2); 
    fit_strength_all->FixParameter(5,sigma_r2); */
    
    // Fix only the temperature parameter 
    fit_strength_all->FixParameter(6,temp);  

    // Set fit range on the spin flip to ensure a narrow resonance following the data points
    //fit_strength_all->SetParLimits(8,0.2,1.);
     
    // Set fit range to something sensible, necessary for the 
    // upper-limit normalization
    //fit_strength_all->SetParLimits(0,16.,18.);
    //fit_strength_all->SetParLimits(3,13.,15.);
    //fit_strength_all->SetParLimits(7,4.5,8.); 
    //fit_strength_all->SetParLimits(8,0.1,1.4); 
    
    cout << endl;
    cout << " Fit of GDR + OCL data, with M1 spin-flip:" << endl;
    OCL_and_GDR_data->Fit(fit_strength_all,"RM+");
    
    fit_strength_all->SetLineColor(kAzure+7);
    fit_strength_all->SetLineWidth(3);
    fit_strength_all->SetLineStyle(1);
    fit_strength_all->Draw("L same");
    
    // Now we have the last parameters in place, we get them and we draw everything nicely.
    E_r1 = fit_strength_all->GetParameter(0);
    Gamma_r1 = fit_strength_all->GetParameter(1);
    sigma_r1 = fit_strength_all->GetParameter(2);
    E_r2 = fit_strength_all->GetParameter(3);
    Gamma_r2 = fit_strength_all->GetParameter(4);
    sigma_r2 = fit_strength_all->GetParameter(5);
    temp = fit_strength_all->GetParameter(6);    
    E_pyg1 = fit_strength_all->GetParameter(7);
    Gamma_pyg1 = fit_strength_all->GetParameter(8);
    sigma_pyg1 = fit_strength_all->GetParameter(9);
    constup1 = fit_strength_all->GetParameter(10);
    constup2 = fit_strength_all->GetParameter(11);
        
    
    TF1 *fit_strength_M1 = new TF1("fit_strength_M1",FitFunctionM1,0.1,30.,5); // SLO M1 + upbend (5 parameters)
    double par3_array[5] = {E_pyg1,Gamma_pyg1,sigma_pyg1,constup1,constup2};
    fit_strength_M1->SetParameters(par3_array);
    fit_strength_M1->FixParameter(0,E_pyg1); 
    fit_strength_M1->FixParameter(1,Gamma_pyg1); 
    fit_strength_M1->FixParameter(2,sigma_pyg1); 
    fit_strength_M1->FixParameter(3,constup1); 
    fit_strength_M1->FixParameter(4,constup2); 
    
    cout << " For plotting only, M1 part:" << endl;
    
    OCL_and_GDR_data->Fit(fit_strength_M1,"RM+");
    
    fit_strength_M1->SetLineColor(kBlack);
    fit_strength_M1->SetLineWidth(1);
    fit_strength_M1->SetLineStyle(3);
    fit_strength_M1->Draw("L same");
    
    // We also plot the E1 component
    TF1 *fit_strength_E1_plot = new TF1("fit_strength_E1_plot",FitFunctionE1,0.1,30.,7); // GDR (7 parameters)
    //Vector for 7 start parameters (GDR)
    double parE1_array[7] = {E_r1,Gamma_r1,sigma_r1,E_r2,Gamma_r2,sigma_r2,temp};
    //                      p0    p1    p2     p3       p4     p5     p6  
    fit_strength_E1_plot->SetParameters(parE1_array);
    
    // Fix all parameters - this is for plotting
    fit_strength_E1_plot->FixParameter(0,E_r1);     
    fit_strength_E1_plot->FixParameter(1,Gamma_r1); 
    fit_strength_E1_plot->FixParameter(2,sigma_r1); 
    fit_strength_E1_plot->FixParameter(3,E_r2); 
    fit_strength_E1_plot->FixParameter(4,Gamma_r2); 
    fit_strength_E1_plot->FixParameter(5,sigma_r2); 
    fit_strength_E1_plot->FixParameter(6,temp);  

    cout << " For plotting only, E1 part:" << endl;    
    OCL_and_GDR_data->Fit(fit_strength_E1_plot,"RM+");
    
    fit_strength_E1_plot->SetLineColor(kBlack);
    fit_strength_E1_plot->SetLineWidth(1);
    fit_strength_E1_plot->SetLineStyle(1);
    fit_strength_E1_plot->Draw("L same");
    cout << " Value of E1 function at 1 MeV: " << fit_strength_E1_plot->Eval(1.0) << endl;

    t.DrawLatex(1.,1.2E-06,"(b)");
    
    
    TLegend *leg = new TLegend(0.2,0.75,0.45,0.93);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.037);
    leg->SetFillColor(0);
    leg->AddEntry(fit_strength_E1_plot," fit E1 ","L");
    leg->AddEntry(fit_strength_M1," fit M1 ","L");
    leg->AddEntry(fit_strength_all," fit sum E1+M1 ","L");
    leg->Draw();
    

    
    
    c1->Update();
    c1->Print("plotstrengthdata_theory_allY.pdf");
    c1->Print("plotstrengthdata_theory_allY.eps");
    c1->Print("plotstrengthdata_theory_allY.png");
    
    // Print E1 and M1 strengths to file, to use in the TALYS calculations
    // REMEMBER that the TALYS functions are given in mb/MeV (Goriely's tables)
    // SO we must convert with the factor const double factor   = 8.6737E-08;	// const. factor in mb^(-1) MeV^(-2)
    FILE *E1file, *M1file;  
    
    cout << endl;
    //cout << " Numbers for the E1 file, in 1/MeV^3: " << endl;
    //E1file = fopen("E1_gsf_90Y_low.txt","w");
    E1file = fopen("E1_gsf_90Y_middle.txt","w");
    //E1file = fopen("E1_gsf_90Y_high.txt","w");
    fprintf(E1file," Z=  39 A=  90\n");
    fprintf(E1file,"  U[MeV]  fM1[mb/MeV]\n");
    double dummy = 0.1;
    for(i=0;i<300;i++){
        fprintf(E1file,"%9.3f%12.3E\n",dummy,(fit_strength_E1_plot->Eval(dummy))/factor);
        //cout << dummy << " " << (fit_strength_E1_plot->Eval(dummy)) << endl;
        dummy += 0.1;
    }
    fclose(E1file);

    //M1file = fopen("M1_gsf_90Y_low.txt","w");
    M1file = fopen("M1_gsf_90Y_middle.txt","w");
    //M1file = fopen("M1_gsf_90Y_high.txt","w");
    fprintf(M1file," Z=  39 A=  90\n");
    fprintf(M1file,"  U[MeV]  fM1[mb/MeV]\n");
    dummy = 0.1;
    for(i=0;i<300;i++){
        fprintf(M1file,"%9.3f%12.3E\n",dummy,(fit_strength_M1->Eval(dummy))/factor);
        //cout << dummy << " " << (fit_strength_M1->Eval(dummy)) << endl;
        dummy += 0.1;
    }
    fclose(M1file);
    
    cout << " Modeled strength functions (rec. norm.) written in TALYS format " << endl;
    //cout << " Modeled strength functions (low norm.) written in TALYS format " << endl;
    //cout << " Modeled strength functions (high norm.) written in TALYS format " << endl;
 
    
    
    
    
    //cout << 11.54E-09*1.8*3400<< endl;
//(pow(0.005,3.)) 
    

    
    
    

}
