#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#define nfiles 100
using namespace std;

//Parameters
//const int nfiles = 100;
const int field_on = 20;
const double vval = 0.1;
const double e_scale = 23.0606;
const double efield = vval * e_scale;
const int nbins = 101;
const int tprobes[] = {field_on, 140, 160, 180, 200, 220, 240, 260, 280, 300};
const int nts = sizeof(tprobes) / sizeof(tprobes[0]);
const int noxy = 512;
const double ocharge = -0.83;
const double ibond_length = 1.0 / 0.58588;

double initiate(double xsize, double x[]){

	double delta = xsize / (nbins - 1);
	for (int i = 0; i<nbins; i++){
		x[i] = i * delta;
	}
	return delta;
}

void get_csts(int &a, int &b, double x[], double y[], double z[], double &xs, double &ys, double &zs){
	
	int nlines = 0;
	string line;
	ifstream indat("fieldoff-1");

	while (getline(indat,line)){
		if (nlines == 3){
			istringstream iss(line);
			iss >> b;
		}
		else if (nlines == 5){
			istringstream iss(line);
			iss >> x[0] >> x[1];
		}
		else if (nlines == 6){
			istringstream iss(line);
			iss >> y[0] >> y[1];
		}
		else if (nlines == 7){
			istringstream iss(line);
			iss >> z[0] >> z[1];
		}
		++nlines;
	}

	a = nlines / (b + 9) - 1;

	x[0] = 0.0;
	xs = (x[1] - x[0]) * 1.000001;
	ys = (y[1] - y[0]) * 1.000001;
	zs = (z[1] - z[0]) * 1.000001;

	indat.close();
}

void read_in(double oxygens[nts][noxy][4], double hydrogens[nts][noxy][2][3], string onoff, int ifile, int nions, double ymin, double zmin, double xsize){

	double xtmp = 0.0, ytmp = 0.0, ztmp = 0.0, utmp = 0.0, q = 0.0;
	int iden = 0, label = 0;
	int i, j, jtmp, olab, hlab;

	stringstream sfile;
	sfile << ifile;
	string filename = "field" + onoff + "-" + sfile.str();
	string line;

	ifstream indat(filename.c_str());

	for (j = 0; j < nts; j++){
		if (j == 0){
			jtmp = tprobes[j] * (nions + 9) + 9;
		}
		else {
			jtmp = (tprobes[j] - tprobes[j-1] - 1) * (nions + 9) + 9;
		}
		for (i = 0; i < jtmp; i++){
			getline(indat, line);
		}
		for (i = 0; i < nions; i++){
			getline(indat, line);
			istringstream iss(line);
			iss >> iden >> label >> q >> xtmp >> ytmp >> ztmp >> utmp; 
			olab = (iden - 1) / 3;
			if (label == 1){
				hlab = (iden + 1)%3;
				hydrogens[j][olab][hlab][0] = xtmp;
				hydrogens[j][olab][hlab][1] = ytmp - ymin;
				hydrogens[j][olab][hlab][2] = ztmp - zmin;
			}
			else if (label == 2){
				oxygens[j][olab][0] = xtmp;
				oxygens[j][olab][1] = ytmp - ymin;
				oxygens[j][olab][2] = ztmp - zmin;
				oxygens[j][olab][3] = utmp;	
				if (onoff == "on"){
					oxygens[j][olab][3] += ocharge*0.5*(xsize-2.0*oxygens[j][olab][0])*efield;
				}
			}
		}
	}

	indat.close();

}

void hist_bin(double hydrogens[nts][noxy][2][3], double oxygens[nts][noxy][4], double xsize, double ysize, double zsize, int onoff, double bulk_pol[2][nbins], double bulk_rho[2][nbins], double bulk_mad[2][nbins], double bulk_mad0[nbins], double bulk_rho0[nbins]){

	int j, k, ii;
	double hbond[3] = {};
	double hbondx, htmp;
	int xlabel;

	for (j = 0; j < nts; j++){
		for (k = 0; k < noxy; k++){

			for (ii = 0; ii < 2; ii++){
				htmp = hydrogens[j][k][ii][1] - oxygens[j][k][1];
				hydrogens[j][k][ii][1] -= int(htmp * 2.0 / ysize);
				htmp = hydrogens[j][k][ii][2] - oxygens[j][k][2];
				hydrogens[j][k][ii][2] -= int(htmp * 2.0 / zsize);
			}

			for (ii = 0; ii < 3; ii++){
				hbond[ii] = 0.5*(hydrogens[j][k][0][ii] + hydrogens[j][k][1][ii]) - oxygens[j][k][ii];
			}
			hbondx = hbond[0] * ibond_length;
			if (hbondx > 1.1){
				cout << "H bond too long!" << endl;
			}

			xlabel = int(oxygens[j][k][0] * nbins / xsize);
			if (j == 0){
				if (onoff == 0){
					bulk_mad0[xlabel] -= oxygens[j][k][3];
				}
				else {
					bulk_mad0[xlabel] += oxygens[j][k][3];
				}
				bulk_rho0[xlabel] += 0.5;
			}
			else {
				bulk_pol[onoff][xlabel] += hbondx;
				bulk_rho[onoff][xlabel] += 1.0;
				bulk_mad[onoff][xlabel] += oxygens[j][k][3];
			}
		}
	}

}

void normalise(double bulk_pol[2][nbins], double bulk_mad[2][nbins], double bulk_rho[2][nbins], double bulk_mad0[nbins], double bulk_rho0[nbins]){

	int i, j;

	for (j = 0; j < nbins; j++){
		for (i = 0; i < 2; i++){
			bulk_pol[i][j] = bulk_pol[i][j] / (bulk_rho[i][j] + int(bulk_rho[i][j] == 0.0));
			bulk_mad[i][j] = bulk_mad[i][j] / (bulk_rho[i][j] + int(bulk_rho[i][j] == 0.0));
			bulk_rho[i][j] = bulk_rho[i][j] / (nfiles * noxy);
		}
		bulk_mad0[j] = bulk_mad0[j] / (bulk_rho0[j] + (bulk_rho0[j] == 0));
	}

}

int main(){

	int nions, nt, ifile;
	double xlims[2] = {}, ylims[2] = {}, zlims[2] = {};
	double xsize, ysize, zsize;
	double xarray[nbins] = {}, deltax;
	double bulk_pol[2][nbins] = {}, bulk_rho[2][nbins] = {}, bulk_mad[2][nbins] = {};
	double bulk_mad0[nbins] = {}, bulk_rho0[nbins] = {};

	double oxygens[nts][noxy][4] = {}, hydrogens[nts][noxy][2][3] = {};

	get_csts(nt, nions, xlims, ylims, zlims, xsize, ysize, zsize);
	deltax = initiate(xlims[1], xarray);

	for (ifile = 0; ifile < nfiles; ifile++){
		cout << ifile+1 << endl;
		read_in(oxygens, hydrogens, "off", ifile+1, nions, ylims[0], zlims[0], xsize);
		hist_bin(hydrogens, oxygens, xsize, ysize, zsize, 0, bulk_pol, bulk_rho, bulk_mad, bulk_mad0, bulk_rho0);

		read_in(oxygens, hydrogens, "on", ifile+1, nions, ylims[0], zlims[0], xsize);
		hist_bin(hydrogens, oxygens, xsize, ysize, zsize, 1, bulk_pol, bulk_rho, bulk_mad, bulk_mad0, bulk_rho0);
	}

	normalise(bulk_pol, bulk_mad, bulk_rho, bulk_mad0, bulk_rho0);

	return 0;
}

