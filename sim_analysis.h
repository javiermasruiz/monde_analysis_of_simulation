#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAttPad.h>
#include "TSystem.h"
#include "TStopwatch.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cstring>

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */

using namespace std;

/*
Autor: Javier Mas Ruiz
Notas: Programa para analizar los datos de salida de la simulacion del 
        detector de neutrones MONDE con geometría cilíndrica.
*/

// Global variables
// Editable 
TString sim_data_input = "sim_data/monde2_out.root"; // simulation data input 
TString data_output_name = "anger_in.txt";  // data input for anger algoritm
TString data_output_anger = "anger_out.txt";  // data output after applying anger algoritm
TString analysis_out = "analysis_out.root";  // analysis output


void readroot_to_txt()
{
    cout << BLUE << " *** *** *** --> BEGIN readroot_to_txt <-- *** *** ***" << RESET << endl;
    cout <<"- Data input file -> " << sim_data_input << endl;
    cout <<"- Data output file-> " << data_output_name << endl;
    cout << " " << endl;

    TFile *finput = new TFile(sim_data_input,"read");
    TDirectory *dir=finput->GetDirectory("ntuple;1");
    dir->ls(); 

    TTree *tree = (TTree*)dir->Get("monde2_ntuple");

    Double_t PosX, PosY;
    Int_t sumPMT1, sumPMT2, sumPMT3, sumPMT4, sumPMT5, sumPMT6;
    Int_t sumPMT7, sumPMT8, sumPMT9, sumPMT10, sumPMT11;
    Int_t sumPMT12, sumPMT13, sumPMT14, sumPMT15, sumPMT16;

    tree->SetBranchAddress("PosX", &PosX);
    tree->SetBranchAddress("PosY", &PosY);
    tree->SetBranchAddress("sumPMT1", &sumPMT1);
    tree->SetBranchAddress("sumPMT2", &sumPMT2);
    tree->SetBranchAddress("sumPMT3", &sumPMT3);
    tree->SetBranchAddress("sumPMT4", &sumPMT4);
    tree->SetBranchAddress("sumPMT5", &sumPMT5);
    tree->SetBranchAddress("sumPMT6", &sumPMT6);
    tree->SetBranchAddress("sumPMT7", &sumPMT7);
    tree->SetBranchAddress("sumPMT8", &sumPMT8);
    tree->SetBranchAddress("sumPMT9", &sumPMT9);
    tree->SetBranchAddress("sumPMT10", &sumPMT10);
    tree->SetBranchAddress("sumPMT11", &sumPMT11);
    tree->SetBranchAddress("sumPMT12", &sumPMT12);
    tree->SetBranchAddress("sumPMT13", &sumPMT13);
    tree->SetBranchAddress("sumPMT14", &sumPMT14);
    tree->SetBranchAddress("sumPMT15", &sumPMT15);
    tree->SetBranchAddress("sumPMT16", &sumPMT16);
    
    ofstream out;
    out.open(data_output_name,ios::out);

    if(out.fail()){
      cout <<" No se pudo abrir el archivo de Salida ";
      exit(1);
    }

    int entries = tree->GetEntries();

    for(int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);
        /*
        cout <<fixed<<showpoint<<setprecision(2)<< PosX << " " << PosY;
        cout << " " << sumPMT1 << " " << sumPMT2 << " " << sumPMT3 << " " << sumPMT4; 
        cout << " " << sumPMT5 << " " << sumPMT6 << " " << sumPMT7 << " " << sumPMT8;
        cout << " " << sumPMT9 << " " << sumPMT10 << " " << sumPMT11 << " " << sumPMT12;
        cout << " " << sumPMT13 << " " << sumPMT14 << " " << sumPMT15 << " " << sumPMT16 << endl;
        */

        out <<fixed<<showpoint<<setprecision(2)<< PosX << " " << PosY;
        out << " " << sumPMT1 << " " << sumPMT2 << " " << sumPMT3 << " " << sumPMT4; 
        out << " " << sumPMT5 << " " << sumPMT6 << " " << sumPMT7 << " " << sumPMT8;
        out << " " << sumPMT9 << " " << sumPMT10 << " " << sumPMT11 << " " << sumPMT12;
        out << " " << sumPMT13 << " " << sumPMT14 << " " << sumPMT15 << " " << sumPMT16 << endl;
    }

    out.close();

    cout << " " << endl;
    cout << "- "<< data_output_name << " --> " << entries << " <-- simulated events" << endl;
    cout << " " << endl;
    // Compile con el flag "-b" si no desea ver la interfaz grafica, ej: root -b macro.C 
    TCanvas *c1 = new TCanvas("c1","Simulated events",500,500);
    tree->Draw("PosY:PosX","","colz");
    c1->SaveAs("simulation.png");
    cout << RED << "*** *** *** --> END <-- *** *** ***" << RESET << endl;
}


void anger(){
    cout << " " << endl;
    cout << BLUE << "*** *** *** --> BEGIN Anger Algorithm <-- *** *** ***" << RESET << endl;

    Int_t threshold = 50;  // Modificar el numero minimo de cuentas a tener en cuenta.
    Int_t good_number_to_run = 4;   // Editar el numero minimo de PMTs con los que se correrá el codigo.
    
    Double_t x, y, X, Y, X0, Y0, dX, dY;
    Int_t pmt_sum, fail_pmt, good_PMT;
    Int_t pmt1, pmt2, pmt3, pmt4, pmt5, pmt6, pmt7, pmt8, pmt9, pmt10, pmt11, pmt12, pmt13, pmt14, pmt15, pmt16;
    Int_t PMT1, PMT2, PMT3, PMT4, PMT5, PMT6, PMT7, PMT8, PMT9, PMT10, PMT11, PMT12, PMT13, PMT14, PMT15, PMT16;

    Double_t x1pos = 30.895, x2pos = 26.191, x3pos = 17.500, x4pos = 6.145, x5pos = -6.145;
    Double_t x6pos = -17.500, x7pos = -26.191, x8pos = -30.891, x9pos = -30.895, x10pos = -26.191;
    Double_t x11pos = -17.500, x12pos = -6.145, x13pos = 6.145, x14pos = 17.500, x15pos = 26.191, x16pos = 30.895;

    Double_t y1pos = 6.145, y2pos = 17.500, y3pos = 26.191, y4pos = 30.895, y5pos = 30.895;
    Double_t y6pos = 26.191, y7pos = 17.500, y8pos = 6.145, y9pos = -6.145, y10pos = -17.500;
    Double_t y11pos = -26.191, y12pos = -30.895, y13pos = -30.895, y14pos = -26.191, y15pos = -17.500, y16pos = -6.145;

    ifstream in;
    ofstream out;
    
    in.open(data_output_name, ios::in); // anger_input_data .txt format
    out.open(data_output_anger, ios::out); // anger_output .txt format

    if(in.fail()){ cout <<" No se pudo abrir el archivo de Entrada "; exit(1); }
    if(out.fail()){ cout <<" No se pudo abrir el archivo de Salida"; exit(1); }

    Int_t entries = 0;
    while (in.good()){  
        in >> x >> y >> pmt1 >> pmt2 >> pmt3 >> pmt4 >> pmt5
           >> pmt6 >> pmt7 >> pmt8 >> pmt9 >> pmt10
           >> pmt11 >> pmt12 >> pmt13 >> pmt14 >> pmt15 >> pmt16;

        if(in.eof()) break;

        X = x;
        Y = y;
        PMT1 = pmt1;
        PMT2 = pmt2;
        PMT3 = pmt3;
        PMT4 = pmt4;
        PMT5 = pmt5;
        PMT6 = pmt6;
        PMT7 = pmt7;
        PMT8 = pmt8;
        PMT9 = pmt9;
        PMT10 = pmt10;
        PMT11 = pmt11;
        PMT12 = pmt12;
        PMT13 = pmt13;
        PMT14 = pmt14;
        PMT15 = pmt15;
        PMT16 = pmt16;
        
        fail_pmt = 0;
        
        if( PMT1 < threshold ){ PMT1 = 0; fail_pmt += 1; }
        if( PMT2 < threshold ){ PMT2 = 0; fail_pmt += 1; }
        if( PMT3 < threshold ){ PMT3 = 0; fail_pmt += 1; }
        if( PMT4 < threshold ){ PMT4 = 0; fail_pmt += 1; }
        if( PMT5 < threshold ){ PMT5 = 0; fail_pmt += 1; }
        if( PMT6 < threshold ){ PMT6 = 0; fail_pmt += 1; }
        if( PMT7 < threshold ){ PMT7 = 0; fail_pmt += 1; }
        if( PMT8 < threshold ){ PMT8 = 0; fail_pmt += 1; }
        if( PMT9 < threshold ){ PMT9 = 0; fail_pmt += 1; }
        if( PMT10 < threshold ){ PMT10 = 0; fail_pmt += 1; }
        if( PMT11 < threshold ){ PMT11 = 0; fail_pmt += 1; }
        if( PMT12 < threshold ){ PMT12 = 0; fail_pmt += 1; }
        if( PMT13 < threshold ){ PMT13 = 0; fail_pmt += 1; }
        if( PMT14 < threshold ){ PMT14 = 0; fail_pmt += 1; }
        if( PMT15 < threshold ){ PMT15 = 0; fail_pmt += 1; }
        if( PMT16 < threshold ){ PMT16 = 0; fail_pmt += 1; } 

        good_PMT = 16 - fail_pmt;

        pmt_sum = PMT1+PMT2+PMT3+PMT4+PMT5+PMT6+PMT7+PMT8+PMT9 + PMT10 + PMT11 + PMT12 + PMT13 + PMT14 + PMT15 + PMT16;
        if(pmt_sum == 0) continue;

        if(good_PMT >= good_number_to_run){

          X0 = (x1pos*PMT1 + x2pos*PMT2 + x3pos*PMT3 + x4pos*PMT4 + x5pos*PMT5 + x6pos*PMT6
                 + x7pos*PMT7 + x8pos*PMT8 + x9pos*PMT9 + x10pos*PMT10 + x11pos*PMT11 
                 + x12pos*PMT12 + x13pos*PMT13 + x14pos*PMT14 + x15pos*PMT15 + x16pos*PMT16)/(pmt_sum*1.0);

          Y0 = (y1pos*PMT1 + y2pos*PMT2 + y3pos*PMT3 + y4pos*PMT4 + y5pos*PMT5 + y6pos*PMT6
                 + y7pos*PMT7 + y8pos*PMT8 + y9pos*PMT9 + y10pos*PMT10 + y11pos*PMT11 
                 + y12pos*PMT12 + y13pos*PMT13 + y14pos*PMT14 + y15pos*PMT15 + y16pos*PMT16)/(pmt_sum*1.0);
        
          dX = fabs(X-X0);
          dY = fabs(Y-Y0);
        }

        
        if(in.good()){
        /*
           cout <<fixed<<showpoint<<setprecision(2)<< X << " " << Y << " " << X0 << " " << Y0 << " " <<  dX  << " " <<  dY;
           cout << " " << PMT1 << " " << PMT2 << " " << PMT3 << " " << PMT4 << " " << PMT5 << " " << PMT6 << " " << PMT7; 
           cout << " " << PMT8 << " " << PMT9 << " " << PMT10 << " " << PMT11 << " " << PMT12 << " " << PMT13 << " " << PMT14;
           cout << " " << PMT15 << " " << PMT16 << "  " << pmt_sum << endl;
        */
           out <<fixed<<showpoint<<setprecision(2)<< X << " " << Y << " " << X0 << " " << Y0 << " " <<  dX  << " " <<  dY;
           out << " " << PMT1 << " " << PMT2 << " " << PMT3 << " " << PMT4 << " " << PMT5 << " " << PMT6 << " " << PMT7; 
           out << " " << PMT8 << " " << PMT9 << " " << PMT10 << " " << PMT11 << " " << PMT12 << " " << PMT13 << " " << PMT14;
           out << " " << PMT15 << " " << PMT16 << "  " << pmt_sum << endl;

          entries++;
          //cout << entries << " -- " << good_PMT  << endl;
        }

    }
    cout << "- "<< data_output_anger << " write -> " << entries << " <- lines." << endl;
    cout << "- Output file --> "<< data_output_anger << endl;

    in.close();
    out.close();
    cout << RED << " *** *** *** --> END <-- *** *** ***" << RESET << endl;
    cout << " " << endl;
}


void plot_anger_root()
{
    cout << BLUE << "*** *** *** --> BEGIN plot_anger_root <-- *** *** ***" << RESET << endl;

    fstream in;
    in.open(data_output_anger, ios::in);

    Double_t x, y, x0, y0, dx, dy;
    Int_t pmt1, pmt2, pmt3, pmt4, pmt5, pmt6;
    Int_t pmt7, pmt8, pmt9, pmt10, pmt11, pmt12;
    Int_t pmt13, pmt14, pmt15, pmt16;
    Int_t pmt_sum;

    TFile *output = new TFile(analysis_out,"recreate");

    TTree *tree = new TTree("tree","tree");

    tree->Branch("x", &x, "x/D");
    tree->Branch("y", &y, "y/D");
    tree->Branch("x0", &x0, "x0/D");
    tree->Branch("y0", &y0, "y0/D");
    tree->Branch("dx", &dx, "dx/D");
    tree->Branch("dy", &dy, "dy/D");

    Int_t entries = 0;
    while(1)
    { 
        in >> x >> y >> x0 >> y0 >> dx >> dy
         >> pmt1 >> pmt2 >> pmt3 >> pmt4
         >> pmt5 >> pmt6 >> pmt7 >> pmt8 >> pmt9
         >> pmt10 >> pmt11 >> pmt12 >> pmt13 >> pmt14
         >> pmt15 >> pmt16 >> pmt_sum;

        if(in.eof()) break;

        //cout << x << " " << y << endl;  

        tree->Fill();

        entries++;
        
    }
    
    cout << "- Analized --> "<< entries << " <-- lines" << endl;


    gStyle->SetPalette(1,0);
    gStyle->SetNumberContours(100);
    gErrorIgnoreLevel = kError;

    TCanvas *c1 = new TCanvas("c1", "Position & Errors",800,800);
    gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
//    gStyle->SetPadTopMargin(0.1);
//    gStyle->SetPadBottomMargin(0.1);

    c1->Divide(2,2,0.01,0.01);

    c1->cd(1);
    //gPad->SetTickx(2);
    tree->Draw("dy:dx","","colz");

    c1->cd(2);
    //gPad->SetTickx(2);
    tree->Draw("dy:dx","","surf3z");

    c1->cd(3);
    tree->Draw("y:x","","colz");

    c1->cd(4);
    tree->Draw("y0:x0","","colz");
    
    c1->SaveAs("PosX_vs_PosY.png");

    tree->Write();

    //output->Write();

    //output->Close();

    in.close();

    cout << "- Output file --> "<< analysis_out << endl;
    cout << RED << " *** *** *** --> END <-- *** *** ***" << RESET << endl;
    cout << " " << endl;

}
