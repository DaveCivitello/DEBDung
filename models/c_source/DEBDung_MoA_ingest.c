/* file DEBDung_MoA_ingest.c */
#include <R.h>
#include <math.h>
static double parms[34];

#define iM parms[0]
#define k parms[1]
#define M parms[2]
#define EM parms[3]
#define Fh parms[4]
#define muD parms[5]
#define DR parms[6]
#define fe parms[7]
#define yRP parms[8]
#define ph parms[9]
#define yPE parms[10]
#define iPM parms[11]
#define eh parms[12]
#define mP parms[13]
#define alpha parms[14]
#define yEF parms[15]
#define LM parms[16]
#define kR parms[17]
#define d0 parms[18]
#define kk parms[19]
#define hb parms[20]
#define theta parms[21]
#define mR parms[22]
#define yVE parms[23]
#define yEF2 parms[24]
#define yEF3 parms[25]
#define yED parms[26]
#define rho parms[27]
#define kR2 parms [28]
#define mR2 parms [29]
#define kk2 parms [30]
#define d02 parms [31]
#define kkM parms [32]
#define d0M parms [33]

/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=34;
odeparms(&N, parms);
}

/* Derivatives and 2 output variables */
void derivs (int *neq, double *t, double *y, double *ydot,
double *yout, int *ip)
{
if (ip[0] <1) error("nout should be at least 1");

/* Dung abundance, in mg */
double TankVol = 0.160;
double DungMass = 1000*y[10]*TankVol; 

/* Toxicant effect. s = stress coefficient */
double s = kkM*(fmax(y[8] - d0M, 0));

/* Other Schito-DEB models use max ingestion rate, half-saturation constant, but here we need max ingestion rate attack rate, a, to set up SU for ingestion */
double a = iM/Fh;

/* Maximum assimilation rate depends on max ingestion rate, iM, and the highest yield of any food type */
double maxYield = fmax(yEF, fmax(yEF2, yEF3));
double aM = maxYield*iM;

/* scaled ingestion rate, fH, which appears in reserve dynamics, needs to account for differences in food yields */
/* So now, this basically needs to be "scaled assimilation rate", since yields decouple these things */
/* Also, this is derived from a substitutable, interacting, sequential supply SU, similar to Lavaud et al. 2014, J Sea Research */
double thetaF = a*y[0]*(a*(DungMass + y[0])*rho + iM)/((a*(DungMass + y[0]) + iM)*(a*y[0]*rho + iM));
double thetaD = a*DungMass*iM/((a*(DungMass + y[0]) + iM)*(a*y[0]*rho + iM));
double fH = (yEF3*thetaF + yED*thetaD)/maxYield/(1 + s);

double Chi = M/(1 + EM);
double fP = y[2]/(y[2] + eh);

double L = y[1];
double LG = fmax(y[1], yout[1]);

double GVOL = pow(LG, 3);
double VOL = pow(L,3);
double SA = pow(L,2);

double Dens = y[5]/(Chi*GVOL);
double kstar = fmin(k + y[5]*alpha, 1);
double g = 1/(yVE*kstar*EM);

double mV = aM*k/(LM*Chi);
double mD = muD*mV;
double rp = Dens*Dens/(ph*ph + Dens*Dens);
double Jec = y[2]*g/(g + y[2])*(aM*SA + yVE*EM*(mV+mR*EM*y[7]+mR2*EM*y[8])*Chi*VOL);

ydot[0] = -iM*thetaF/(1 + s)*SA;
ydot[10] = -iM*thetaD/(1 + s)*SA/1000/TankVol;
ydot[1] = yVE/(3*Chi*SA)*(kstar*Jec - (mV+mR*EM*y[7]+mR2*EM*y[8])*Chi*VOL);
ydot[2] = aM/(Chi*EM*L)*(fH - y[2]) - iPM*y[5]*fP/(EM*Chi*VOL);
ydot[5] = yPE*iPM*fP*(1 - rp)*y[5] - mP*y[5];
ydot[6] = fmax(yRP*yPE*iPM*fP*rp*y[5],0);
if(y[3] < DR){
  ydot[3] = (1 - kstar)*Jec - mD*y[3];
  ydot[4] = 0;}else{
  ydot[3] = fmin(0, (1 - kstar)*Jec - mD*DR);
  ydot[4] = fmax((1 - kstar)*Jec - mD*DR, 0);}
ydot[7] = theta/(Chi*VOL)*ydot[6] + kR*fmax(1-y[2],0) - kR*y[7] - 3*y[7]*ydot[1]/L;
ydot[8] = kR2*(y[10] - y[8]) - 3*y[8]*ydot[1]/L;
ydot[9] = kk*fmax(y[7] - d0, 0) + kk2*fmax(y[8] - d02, 0) + hb;
if(y[2] <= 0){
  ydot[0] = 0;
  ydot[1] = 0;
  ydot[2] = 0;
  ydot[3] = 0;
  ydot[4] = 0;
  ydot[5] = 0;
  ydot[6] = 0;
  ydot[7] = 0;
  ydot[9] = hb;
  ydot[10] = 0;}

  yout[0] = exp(-y[9]);
  yout[1] = LG;

}

/* END file DEBDung_MoA_ingest.c */ 
