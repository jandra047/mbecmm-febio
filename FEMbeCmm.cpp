#include "stdafx.h"
#include "FEMbeCmm.h"
#include "FECore/FEAnalysis.h"					// To get end time
#include "FECore/FEModel.h"						// To get current time
#include <iostream>								// To use cin.get()
#include <sstream>
#include <signal.h>
#define _USE_MATH_DEFINES						// To introduce pi constant (1/2)
#include <math.h>								// To introduce pi constant (2/2)
#include "FECore/log.h"							// To print to log file and/or screen

GRMaterialPoint::GRMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}

// Might be unneccessary
FEMaterialPoint* GRMaterialPoint::Copy()
{
	GRMaterialPoint* pt = new GRMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

void GRMaterialPoint::Init()
{
	FEMaterialPoint::Init();

	m_Jo = 1;
	m_svo = 0;
	m_smo.zero();
	m_sco.zero();
	m_Fio.unit();
	m_Jh = 1;
	m_Fih.unit();

	m_phic = 0;
	m_Iemax = 0;
}

void GRMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPoint::Serialize(ar);
	ar & m_Jo & m_svo & m_smo & m_sco & m_Fio & m_Jh & m_Fih & m_phic & m_Iemax;
}

// Added by David Li (Sept 2021) to circumvent Error C4716 FEMbeCMM::Tangent must return a value
tens4ds FEMbeCmm::Tangent(FEMaterialPoint & pt)
{
	return tens4ds();
}

BEGIN_FECORE_CLASS(FEMbeCmm, FEElasticMaterial)
	ADD_PARAMETER(rIo,         FE_RANGE_GREATER(0.0), "rIo");
	ADD_PARAMETER(lo,          FE_RANGE_GREATER(0.0), "lo");
	ADD_PARAMETER(endtime,     FE_RANGE_GREATER(0.0), "endtime");
	ADD_PARAMETER(partialtime, FE_RANGE_GREATER(0.0), "partialtime");
	ADD_PARAMETER(Jdep,        FE_RANGE_GREATER(0.0), "Jdep");
	ADD_PARAMETER(bulkLM,      FE_RANGE_GREATER(0.0), "bulkLM");
	ADD_PARAMETER(imper,       FE_RANGE_GREATER_OR_EQUAL(0.0), "imper");
	ADD_PARAMETER(hwaves,      FE_RANGE_GREATER(0.0), "hwaves");
	ADD_PARAMETER(phieo,       FE_RANGE_CLOSED(0.0, 1.0), "phieo");
	ADD_PARAMETER(phimo,       FE_RANGE_CLOSED(0.0, 1.0), "phimo");
	ADD_PARAMETER(phico,       FE_RANGE_CLOSED(0.0, 1.0), "phico");
	ADD_PARAMETER(betat,       FE_RANGE_CLOSED(0.0, 1.0), "betat");
	ADD_PARAMETER(betaz,       FE_RANGE_CLOSED(0.0, 1.0), "betaz");
	ADD_PARAMETER(alphao,      FE_RANGE_DONT_CARE(), "alphao");
	ADD_PARAMETER(mub,         FE_RANGE_GREATER(0.0), "ce");
	ADD_PARAMETER(cm,          FE_RANGE_GREATER(0.0), "cm1");
	ADD_PARAMETER(dm,          FE_RANGE_GREATER(0.0), "cm2");
	ADD_PARAMETER(ccb,         FE_RANGE_GREATER(0.0), "cc1");
	ADD_PARAMETER(dc,          FE_RANGE_GREATER(0.0), "cc2");
	ADD_PARAMETER(Get,         FE_RANGE_GREATER(0.0), "Get");
	ADD_PARAMETER(Gez,         FE_RANGE_GREATER(0.0), "Gez");
	ADD_PARAMETER(Gm,          FE_RANGE_GREATER(0.0), "Gm");
	ADD_PARAMETER(Gcb,         FE_RANGE_GREATER(0.0), "Gc");
	ADD_PARAMETER(Tmaxb,       FE_RANGE_GREATER_OR_EQUAL(0.0), "Tmax");
	ADD_PARAMETER(CB,          FE_RANGE_GREATER(0.0), "CB");
	ADD_PARAMETER(lamM,        FE_RANGE_GREATER(0.0), "lamM");
	ADD_PARAMETER(lam0,        FE_RANGE_GREATER(0.0), "lam0");
	ADD_PARAMETER(aexp,        FE_RANGE_GREATER_OR_EQUAL(0.0), "aexp");
	ADD_PARAMETER(etab,        FE_RANGE_CLOSED(0.0, 1.0), "eta");
	ADD_PARAMETER(KsKib,       FE_RANGE_DONT_CARE(), "KsKi");
	ADD_PARAMETER(infl,        FE_RANGE_CLOSED(0.0, 1.0), "infl");
	ADD_PARAMETER(KfKicen,     FE_RANGE_DONT_CARE(), "KfKi");
	ADD_PARAMETER(deltab,      FE_RANGE_CLOSED(0.0, 1.0), "delta");
	ADD_PARAMETER(asym,        FE_RANGE_CLOSED(0.0, 1.0), "asym");
	ADD_PARAMETER(zod,         FE_RANGE_GREATER(0.0), "zod");
	ADD_PARAMETER(nuz,         FE_RANGE_GREATER(0.0), "nuz");
	ADD_PARAMETER(tod,         FE_RANGE_GREATER(0.0), "tod");
	ADD_PARAMETER(nut,         FE_RANGE_GREATER(0.0), "nut");
	ADD_PARAMETER(zo,          FE_RANGE_GREATER(0.0), "zo");
	ADD_PARAMETER(to,          FE_RANGE_CLOSED(0.0, 2.0), "to");
	ADD_PARAMETER(insmu,       FE_RANGE_CLOSED(-1.0, 1.0), "insce");
	ADD_PARAMETER(inscc,       FE_RANGE_CLOSED(-1.0, 1.0), "inscc1");
	ADD_PARAMETER(insTmax,     FE_RANGE_CLOSED(-1.0, 1.0), "insTmax");
	ADD_PARAMETER(insdelta,    FE_RANGE_CLOSED(0.0, 1.0), "insdelta");
	ADD_PARAMETER(insGc,       FE_RANGE_CLOSED(-1.0, 1.0), "insGc");
	ADD_PARAMETER(agimu,       FE_RANGE_CLOSED(0.0, 1.0), "agice");
	ADD_PARAMETER(etacen,      FE_RANGE_CLOSED(0.0, 1.0), "etacen");
	ADD_PARAMETER(KsKicen,     FE_RANGE_DONT_CARE(), "KsKicen");
	ADD_PARAMETER(adho,        FE_RANGE_CLOSED(-1.0, 1.0), "adho");
	ADD_PARAMETER(m_insult,    FE_RANGE_CLOSED(0.0, 1.0), "insult");
END_FECORE_CLASS();

void FEMbeCmm::StressTangent(FEMaterialPoint& mp, mat3ds& stress, tens4dmm& tangent)
{
	// The FEMaterialPoint classes are stored in a linked list. The specific material
	// point data needed by this function can be accessed using the ExtractData member.
	// In this case, we want to get FEElasticMaterialPoint data since it stores the deformation
	// information that is needed to evaluate the stress.
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	GRMaterialPoint& pt = *mp.ExtractData<GRMaterialPoint>();
	
	// We'll need the deformation gradient and its determinant in this function.
	// Note that we don't take the determinant of F directly (using mat3d::det)
	// but instead use the m_J member variable of FEElasticMaterialPoint.
	//const mat3d &F = et.m_F;
	//const double J = et.m_J;
	mat3d &F = et.m_F;
	double J = et.m_J;

	const double eps = std::numeric_limits<double>::epsilon();

	// Get current and end times
	const double t = GetFEModel()->GetTime().currentTime;

	//const double endtime = 11.0;							// 11.0 | 31.0-32.0 (TEVG)
	//const double partialtime = endtime;			// partialtime <= endtime | 10.0 | 10.4 (for TI calculation)
	const double sgr = min(t, partialtime);		// G&R time s, min(t,partialtime) | min(t,9.0)

	// Retrieve material position
	const vec3d X = pt.m_r0;
	
	//const double imper = 0.00;					// imper > 0 for TORTUOSITY (see Matlab script <NodesElementsAsy.m>) | 0.00 | 20.0
	//const double rIo = 0.6468;					// 0.6468 | 0.5678
	//const double hwaves = 2.0;
	//const double lo = 15.0;
	const vec3d Xcl = {0.0, imper/100.0*rIo*sin(hwaves*M_PI*X.z/lo), X.z};		// Center line
	
	vec3d NX = {X.x-Xcl.x, X.y-Xcl.y, X.z-Xcl.z};								// Radial vector
	const double ro = sqrt(NX*NX);
	NX /= ro;

	// Retrieve local element basis directions
	vec3d N[3];

	// Point-wise, consistent with mesh generated with Matlab script <NodesElementsAsy.m>
	N[2] = {0.0, imper/100.0*rIo*hwaves*M_PI/lo*cos(hwaves*M_PI*X.z/lo), 1.0}; N[2] = N[2]/sqrt(N[2]*N[2]);		// Axial = d(Xcl)/d(z)
	N[1] = {-NX.y, NX.x, NX.z};																					// Circumferential
	N[0] = N[2]^N[1];

	// Element-wise, from input file (doesn't work, m_Q variable deprecated)
	// N[2] = pt.m_Q.col(0); N[1] = pt.m_Q.col(1); N[0] = pt.m_Q.col(2);							// Axial, circumferential, radial

	//const double phieo = 0.34;								// 0.34 (CMAME | KNOCKOUTS) | 1.00 (TEVG) | 1.0/3.0 (TEVG)
	//const double phimo = 0.5*(1.0 - phieo);
	//const double phico = 0.5*(1.0 - phieo);

	//const double eta = 1.0;									// 1.0 | 1.0/3.0 (for uniform cases) | 0.714

	//const double mu = 89.71;
	//const double Get = 1.90;
	//const double Gez = 1.62;

	//double alpha = 0.522;								// Original orientation of diagonal collagen | 0.522 (CMAME | KNOCKOUTS) | 0.8713 (TEVG)
	double alpha = alphao;

	// Original homeostatic parameters (adaptive)

	// Passive
	//const double cm = 261.4;									// 261.4 (CMAME | KNOCKOUTS) | 46.61 (TEVG)
	//const double dm = 0.24;
	//const double Gm = 1.20;
	//const double cc = 234.9;									// 234.9 (CMAME | KNOCKOUTS) | 328.475 (TEVG)
	//const double dc = 4.08;
	//const double Gc = 1.25;

	// Orientation fractions for collagen
	//const double betat = 0.056;
	//const double betaz = 0.067;
	const double betad = 0.5*(1.0 - betat - betaz);

	vec3d  Np = N[1] * sin(alpha) + N[2] * cos(alpha);		// Original diagonal fiber direction
	vec3d  Nn = N[1] * sin(alpha) - N[2] * cos(alpha);		// idem for symmetric

	// Active
	//const double Tmax = 250.0 * 0.0;							// 250.0 | 50.0 | 150.0 (for uniform cases, except for contractility -> 250)
	//const double lamM = 1.1;
	//const double lam0 = 0.4;
	//const double CB = sqrt(log(2.0));							// Such that (1-exp(-CB^2)) = 0.5
	const double CS = 0.5*CB * 1.0;							// Such that (1-exp( -C^2)) = 0.0 for lt = 1/(1+CB/CS)^(1/3) = 0.7 and (1-exp(-C^2)) = 0.75 for lt = 2.0
	
	//const double KsKi = 0.35;
	const double EPS = 1.0 + (1.0 - 1.0)*(sgr-1.0)/(partialtime-1.0);
	
	//const double KfKi = 1.0;
	//const double inflam = 0.0*(sgr-1.0)/(partialtime-1.0);
	double inflam = infl*(sgr-1.0)/(partialtime-1.0);
	
	//const double aexp = 1.0;									// 1.0 (KNOCKOUTS | TEVG) | 0.0 (CMAME | TORTUOSITY)

	//const double delta = 0.0;

	// Insults occurring when t > 1
	double azi = acos(-NX.y);   // Azimuth wrt axis -Y
	double mu = mub;			// Baseline mu
	double cc = ccb;			// Baseline cc
	double Tmax = Tmaxb;		// Baseline Tmax
	double delta = deltab;		// Baseline mechanosensing
	double Gc = Gcb;			// Baseline collagen deposition stretch
	double eta = etab;			// Baseline eta
	double KsKi = KsKib;		// Baseline KsKi
	double KfKi = 0.0;          // Baseline KfKi
	
	if (t > 1.0 + eps) {

		// mu
		double muout = (1.0 - agimu)*mub;
		double muin = (1.0 - agimu - insmu)*mub;
		mu = muout + (muin - muout)*exp(-pow(abs((X.z-zo)/(lo/2.0/zod)),nuz))*exp(-asym*pow(abs((azi-to*M_PI)/(M_PI/tod)),nut))*(sgr-1.0)/(partialtime-1.0);
		//mu = mub*(1.0 - agimu - insmu*m_insult(mp)*(sgr-1.0)/(partialtime-1.0));

		// cc
		double ccout = 1.0*ccb;
		double ccin = (1.0 - inscc)*ccb;
		cc = ccout + (ccin - ccout)*exp(-pow(abs((X.z-zo)/(lo/2.0/zod)),nuz))*exp(-asym*pow(abs((azi-to*M_PI)/(M_PI/tod)),nut))*(sgr-1.0)/(partialtime-1.0);
		//cc = ccb*(1.0 - inscc*m_insult(mp)*(sgr-1.0)/(partialtime-1.0));

		// Tmax
		double Tmaxout = 1.0*Tmaxb;
		double Tmaxin = (1.0 - insTmax)*Tmaxb;
		Tmax = Tmaxout + (Tmaxin - Tmaxout)*exp(-pow(abs((X.z-zo)/(lo/2.0/zod)),nuz))*(sgr-1.0)/(partialtime-1.0);
		//Tmax = Tmaxb*(1.0 - insTmax*m_insult(mp)*(sgr-1.0)/(partialtime-1.0));

		// delta
		double deltaout = deltab;
		double deltain = deltab + insdelta;
		delta = deltaout + (deltain - deltaout)*exp(-pow(abs((X.z-zo)/(lo/2.0/zod)),nuz))*exp(-asym*pow(abs((azi-to*M_PI)/(M_PI/tod)),nut))*(sgr-1.0)/(partialtime-1.0);
		//delta = deltab + insdelta*m_insult(mp)*(sgr-1.0)/(partialtime-1.0);

		// Gc
		double Gcout = 1.0*Gcb;
		double Gcin = (1.0 - insGc)*Gcb;
		Gc = Gcout + (Gcin - Gcout)*exp(-pow(abs((X.z-zo)/(lo/2.0/zod)),nuz))*exp(-asym*pow(abs((azi-to*M_PI)/(M_PI/tod)),nut))*(sgr-1.0)/(partialtime-1.0);
		//Gc = Gcb * (1.0 - insGc*m_insult(mp)*(sgr-1.0)/(partialtime-1.0));

		// eta
		double etaout = etab;
		double etain = etacen;
		eta = etaout + (etain - etaout)*exp(-pow(abs((X.z-zo)/(lo/2.0/zod)),nuz));
		//eta = etab + (etacen - etab)*m_insult(mp);

		// KsKi
		double KsKiout = KsKib;
		double KsKiin = KsKicen;
		KsKi = KsKiout + (KsKiin - KsKiout)*exp(-pow(abs((X.z-zo)/(lo/2.0/zod)),nuz));
		//KsKi = KsKib + (KsKicen - KsKib)*m_insult(mp);

		// KfKi
		double KfKiout = 0.0;
		double KfKiin = KfKicen;
		KfKi = KfKiout + (KfKiin - KfKiout)*exp(-pow(abs((X.z-zo)/(lo/2.0/zod)),nuz));
	}

	// Compute U from polar decomposition of deformation gradient tensor
	mat3ds U; mat3d R; F.right_polar(R, U);

	double eigenval[3]; vec3d eigenvec[3];
	U.eigen2(eigenval, eigenvec);

	// Right Cauchy-Green tensor and its inverse
	const mat3ds C = et.RightCauchyGreen();
	const mat3ds Ci = C.inverse();

	// Ge from spectral decomposition
	const mat3ds Ge = 1.0/Get/Gez*dyad(N[0])+Get*dyad(N[1])+Gez*dyad(N[2]);
	
	// Stress for elastin
	mat3ds Se = (phieo*mu*Ge*Ge).sym();	// const					// phieo*Ge*Sehat*Ge = phieo*Ge*(mu*I)*Ge
	
	// Computation of the second Piola-Kirchhoff stress
	mat3ds S;

	// Define identity tensor and some useful dyadic products of the identity tensor
	const mat3dd  I(1.0);
	const tens4ds IxI = dyad1s(I);
	const tens4ds IoI = dyad4s(I);
	const tens4dmm IxIss = tens4dmm(IxI);							// IxI in tens4dmm form
	const tens4dmm IoIss = tens4dmm(IoI);							// IoI in tens4dmm form

	// Spatial moduli for elastin
	tens4ds ce(0.0);								// phieo/J*(FcF:GecGe:Cehat:GecGe:FTcFT) = phieo/J*(FcF:GecGe:0:GecGe:FTcFT)

	// Computation of spatial moduli
	tens4dmm css;    // Tangent tensor
	mat3ds sfpro;

	if (t <= 1.0 + eps) {
		// Compute stress
		//const double Jdep = 0.9999;
		//const double lm = 1.0e3*mu;
		const double lm = bulkLM*mu;
		
		const double lt = (F*N[1]).norm();
		const double lz = (F*N[2]).norm();
		const double lp = (F*Np).norm();
		const double ln = (F*Nn).norm();
		
		const double lmt2 = (Gm*lt)*(Gm*lt);
		const double lct2 = (Gc*lt)*(Gc*lt);
		const double lcz2 = (Gc*lz)*(Gc*lz);
		const double lcp2 = (Gc*lp)*(Gc*lp);
		const double lcn2 = (Gc*ln)*(Gc*ln);
		
		//et.m_v.x = sqrt(C.dotdot(dyad(N[0])));					// Radial stretch
		//et.m_v.y = sqrt(C.dotdot(dyad(N[1])));					// Circ. stretch
		//et.m_v.z = sqrt(C.dotdot(dyad(N[2])));					// Axial stretch

		// Write strain energy density and circ, axial stiffness to components of acceleration output for plotting
		// Strain energy density
		double sede = phieo*(mu/2.0*((F*Ge).dotdot(F*Ge)-3.0));
		double sedm = phimo*(cm/(4.0*dm)*(exp(dm*(lmt2-1.0)*(lmt2-1.0))-1.0));
		double sedc = phico*(cc/(4.0*dc)*(exp(dc*(lct2-1.0)*(lct2-1.0))-1.0)*betat+
							 cc/(4.0*dc)*(exp(dc*(lcz2-1.0)*(lcz2-1.0))-1.0)*betaz+
							 cc/(4.0*dc)*(exp(dc*(lcp2-1.0)*(lcp2-1.0))-1.0)*betad+
							 cc/(4.0*dc)*(exp(dc*(lcn2-1.0)*(lcn2-1.0))-1.0)*betad);
		
		et.m_v.y = phimo;
		et.m_v.z = phico;

		et.m_a.x = sede + sedm + sedc;
		
		// Circumferential and axial stiffness components
		et.m_a.y = phieo*(2.0*mu*(((F*Ge).transpose())*(F*Ge)).dotdot(dyad(N[1])))+
				   phimo*(2.0*cm*exp(dm*(lmt2-1.0)*(lmt2-1.0))*((lmt2-1.0)+(1.0+2.0*dm*(lmt2-1.0)*(lmt2-1.0))*pow(Gm,2))*pow(Gm,2))+
				   phico*(2.0*cc*exp(dc*(lct2-1.0)*(lct2-1.0))*((lct2-1.0)+(1.0+2.0*dc*(lct2-1.0)*(lct2-1.0))*pow(Gc,2))*pow(Gc,2)*betat+
						  2.0*cc*exp(dc*(lcp2-1.0)*(lcp2-1.0))*((lcp2-1.0)+(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,2)*pow(sin(alpha),2))*pow(Gc,2)*pow(sin(alpha),2)*betad+
						  2.0*cc*exp(dc*(lcn2-1.0)*(lcn2-1.0))*((lcn2-1.0)+(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,2)*pow(sin(alpha),2))*pow(Gc,2)*pow(sin(alpha),2)*betad);    // c_tttt
		et.m_a.z = phieo*(2.0*mu*(((F*Ge).transpose())*(F*Ge)).dotdot(dyad(N[2])))+
				   phico*(2.0*cc*exp(dc*(lcz2-1.0)*(lcz2-1.0))*((lcz2-1.0)+(1.0+2.0*dc*(lcz2-1.0)*(lcz2-1.0))*pow(Gc,2))*pow(Gc,2)*betaz+
						  2.0*cc*exp(dc*(lcp2-1.0)*(lcp2-1.0))*((lcp2-1.0)+(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,2)*pow(cos(alpha),2))*pow(Gc,2)*pow(cos(alpha),2)*betad+
						  2.0*cc*exp(dc*(lcn2-1.0)*(lcn2-1.0))*((lcn2-1.0)+(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,2)*pow(cos(alpha),2))*pow(Gc,2)*pow(cos(alpha),2)*betad);    // c_zzzz

		// Passive stress
		const mat3ds Sm = (cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm)*dyad(N[1]));
		const mat3ds Sc = (cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc)*dyad(N[1])*betat +
						   cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc)*dyad(N[2])*betaz +
					 	   cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc)*dyad( Np )*betad +
						   cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc)*dyad( Nn )*betad);

		// Active stress
		mat3ds Sa = Tmax*(1.0-exp(-CB*CB))*(1.0-pow((lamM-1.0)/(lamM-lam0),2))*(lt*lt)*dyad(N[1]); // const
		
		mat3ds Sx = Se + phimo * Sm + phico * Sc + phimo * Sa; // const
		
		S = Sx + Ci*lm*log(Jdep*J);
		
		mat3d u(U);
		
		pt.m_Jo    = J;
		pt.m_svo   = 1.0/3.0/J*S.dotdot(C);
		pt.m_smo   = 1.0/J*(u*(Sm*u)).sym();
		pt.m_sco   = 1.0/J*(u*(Sc*u)).sym();
		pt.m_Fio   = F.inverse();
		pt.m_Jh    = pt.m_Jo;
		pt.m_Fih   = pt.m_Fio;
		pt.m_phic  = phico;

		// Compute hyperelastic tangent
		const mat3ds tent = dyad(F*N[1]);
		const mat3ds tenz = dyad(F*N[2]);
		const mat3ds tenp = dyad(F*Np);
		const mat3ds tenn = dyad(F*Nn);

		// Passive
		tens4ds cf = phimo*(2.0*cm*(1.0+2.0*dm*(lmt2-1.0)*(lmt2-1.0))*exp(dm*(lmt2-1.0)*(lmt2-1.0))*pow(Gm,4)*dyad1s(tent))      +
					 phico*(2.0*cc*(1.0+2.0*dc*(lct2-1.0)*(lct2-1.0))*exp(dc*(lct2-1.0)*(lct2-1.0))*pow(Gc,4)*dyad1s(tent)*betat +
							2.0*cc*(1.0+2.0*dc*(lcz2-1.0)*(lcz2-1.0))*exp(dc*(lcz2-1.0)*(lcz2-1.0))*pow(Gc,4)*dyad1s(tenz)*betaz +
							2.0*cc*(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*exp(dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,4)*dyad1s(tenp)*betad +
							2.0*cc*(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*exp(dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,4)*dyad1s(tenn)*betad);

		// Active
		tens4ds ca = phimo*2.0*Tmax*(1.0-exp(-CB*CB))*(1.0-pow((lamM-1.0)/(lamM-lam0),2))*dyad1s(tent);

		cf /= J; ca /= J;

		tens4ds c = ce + cf + ca;

		c += lm/J*(IxI-2.0*log(Jdep*J)*IoI);

		css = tens4dmm(c);		// c in tens4dmm form
	}
	else if (t <= partialtime + eps) {
		// Compute stress
		double    Jo = pt.m_Jo; // const
		double   svo = pt.m_svo; // const
		mat3ds  &smo = pt.m_smo;
		mat3ds  &sco = pt.m_sco;
		mat3d    Fio = pt.m_Fio; // const
		double &phic = pt.m_phic;
		
		phic = phico;																// Initial guess
		double dRdc = J/Jo*(1.0+phimo/phico*eta*pow(J/Jo*phic/phico,eta-1.0));		// Initial tangent d(R)/d(phic)
		double Rphi = phieo+phimo*pow(J/Jo*phic/phico,eta)+J/Jo*phic-J/Jo;			// Initial residue
		do {																		// Local iterations to obtain phic
			phic = phic-Rphi/dRdc;													// phic
			dRdc = J/Jo*(1.0+phimo/phico*eta*pow(J/Jo*phic/phico,eta-1.0));			// Tangent
			Rphi = phieo+phimo*pow(J/Jo*phic/phico,eta)+J/Jo*phic-J/Jo;				// Update residue
		} while (abs(Rphi) > sqrt(eps));											// && abs(Rphi/Rphi0) > sqrt(eps) && j<10
		phic = phic-Rphi/dRdc;														// Converge phase -> phic (updated in material point memory)

		double phim = phimo/(J/Jo)*pow(J/Jo*phic/phico,eta);	// phim from <J*phim/phimo=(J*phic/phico)^eta>
		double phie = phieo/(J/Jo);
		// const const

		// Recompute remodeled original stresses for smc and collagen (from remodeled natural configurations)

		const double lto = (Fio.inverse()*N[1]).norm();
		const double lzo = (Fio.inverse()*N[2]).norm();
		const double lpo = (Fio.inverse()*Np).norm();					// Original referential stretch for deposition stretch calculation
		const double lno = (Fio.inverse()*Nn).norm();					// idem for symmetric

		const double lmt2 = (Gm*lto)*(Gm*lto);
		const double lct2 = (Gc*lto)*(Gc*lto);
		const double lcz2 = (Gc*lzo)*(Gc*lzo);
		const double lcp2 = (Gc*lpo)*(Gc*lpo);						// Deposition stretch calculation (computational purposes)
		const double lcn2 = (Gc*lno)*(Gc*lno);						// idem for symmetric

		const double lr = (F*(Fio*N[0])).norm();						// lr -> 1 for F -> Fo
		const double lt = (F*(Fio*N[1])).norm();						// lt -> 1 for F -> Fo
		const double lz = (F*(Fio*N[2])).norm();						// lz -> 1 for F -> Fo

		alpha = atan(tan(alpha)*pow(lt/lz,aexp));				// Update alpha
		Np = N[1]*sin(alpha)+N[2]*cos(alpha);					// Update diagonal fiber vector
		Nn = N[1]*sin(alpha)-N[2]*cos(alpha);					// idem for symmetric
		
		// Strain energy density
		double sede=phie*(mu/2.0*((F*Ge).dotdot(F*Ge)-3.0));
		double sedm=phim*(cm/(4.0*dm)*(exp(dm*(lmt2-1.0)*(lmt2-1.0))-1.0));
		double sedc=phic*(cc/(4.0*dc)*(exp(dc*(lct2-1.0)*(lct2-1.0))-1.0)*betat+
						  cc/(4.0*dc)*(exp(dc*(lcz2-1.0)*(lcz2-1.0))-1.0)*betaz+
						  cc/(4.0*dc)*(exp(dc*(lcp2-1.0)*(lcp2-1.0))-1.0)*betad+
						  cc/(4.0*dc)*(exp(dc*(lcn2-1.0)*(lcn2-1.0))-1.0)*betad);
		
		et.m_v.y = phim;
		et.m_v.z = phic;

		et.m_a.x = sede + sedm + sedc;
		
		// Circumferential and axial stiffness components
		et.m_a.y = phie*(2.0*mu*(((F*Ge).transpose())*(F*Ge)).dotdot(dyad(N[1])))+
				   phim*(2.0*cm*exp(dm*(lmt2-1.0)*(lmt2-1.0))*((lmt2-1.0)+(1.0+2.0*dm*(lmt2-1.0)*(lmt2-1.0))*pow(Gm,2))*pow(Gm,2))+
				   phic*(2.0*cc*exp(dc*(lct2-1.0)*(lct2-1.0))*((lct2-1.0)+(1.0+2.0*dc*(lct2-1.0)*(lct2-1.0))*pow(Gc,2))*pow(Gc,2)*betat+
						 2.0*cc*exp(dc*(lcp2-1.0)*(lcp2-1.0))*((lcp2-1.0)+(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,2)*pow(sin(alpha),2))*pow(Gc,2)*pow(sin(alpha),2)*betad+
						 2.0*cc*exp(dc*(lcn2-1.0)*(lcn2-1.0))*((lcn2-1.0)+(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,2)*pow(sin(alpha),2))*pow(Gc,2)*pow(sin(alpha),2)*betad);	// c_tttt
		et.m_a.z = phie*(2.0*mu*(((F*Ge).transpose())*(F*Ge)).dotdot(dyad(N[2])))+
				   phic*(2.0*cc*exp(dc*(lcz2-1.0)*(lcz2-1.0))*((lcz2-1.0)+(1.0+2.0*dc*(lcz2-1.0)*(lcz2-1.0))*pow(Gc,2))*pow(Gc,2)*betaz+
						 2.0*cc*exp(dc*(lcp2-1.0)*(lcp2-1.0))*((lcp2-1.0)+(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,2)*pow(cos(alpha),2))*pow(Gc,2)*pow(cos(alpha),2)*betad+
						 2.0*cc*exp(dc*(lcn2-1.0)*(lcn2-1.0))*((lcn2-1.0)+(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,2)*pow(cos(alpha),2))*pow(Gc,2)*pow(cos(alpha),2)*betad);	// c_zzzz

		// Passive stress
		const mat3ds Smo = (cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm)*dyad(N[1]));
		const mat3ds Sco = (cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc)*dyad(N[1])*betat +
							cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc)*dyad(N[2])*betaz +
							cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc)*dyad( Np )*betad +
							cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc)*dyad( Nn )*betad );

		// Active stress
		mat3ds Sao = Tmax*(1.0-exp(-CB*CB))*(1.0-pow((lamM-1.0)/(lamM-lam0),2))*(lto*lto)*dyad(N[1]); // const
		
		mat3ds Uo; mat3d Ro; (Fio.inverse()).right_polar(Ro,Uo);	// Uo from polar decomposition
		mat3d  uo(Uo);
		
		smo = 1.0/Jo*(uo*(Smo*uo)).sym();
		sco = 1.0/Jo*(uo*(Sco*uo)).sym();

		mat3ds sao = 1.0/Jo*(uo*(Sao*uo)).sym(); // const
		
		// Compute current stresses (originally from Stress function)
		
		double rIrIo = ro/rIo*lt-(ro-rIo)/rIo*lr;				// rIrIo -> rIorIo = 1 for F -> Fo

		mat3ds sNm = phim*smo;									// phim*smhato = phim*smo
		mat3ds sNc = phic*sco;									// phic*schato = phic*sco

		const mat3ds sNf = sNm + sNc;
		
		double Cratio = CB-CS*(EPS*pow(rIrIo,-3)-1.0); // const
		mat3ds sNa; sNa.zero();
		if (Cratio>0) sNa = phim*(1.0-exp(-Cratio*Cratio))/(1.0-exp(-CB*CB))*sao;

		mat3ds Ui = U.inverse();  // const          					// Inverse of U
		mat3d  ui(Ui); // const

		// 2nd P-K stresses
		mat3ds Sf = J*(ui*sNf*ui).sym(); // const					// J*Ui*sNf*Ui
		mat3ds Sa = J*(ui*sNa*ui).sym(); // const						// J*Ui*sNa*Ui
		mat3ds Sx = Se + Sf + Sa; // const

		//double p = 1.0/3.0/J*Sx.dotdot(C) - svo/(1.0-delta)*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0)-KfKi*inflam); // const		// Ups = 1 -> p
		double p = 1.0/3.0/J*Sx.dotdot(C) - svo*(1.0 + adho*(sgr-1.0)/(partialtime-1.0))/(1.0-delta)*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0)-KfKi*inflam);		// Ups = 1 -> p
		
		S = Sx - J*p*Ci;
		
		pt.m_Jh  = J;
		pt.m_Fih = F.inverse();

		//// Redundant code from Tangent, some of which is needed, added by David Li (Dec 2021)
		//Jo = pt.m_Jo;
		//svo = pt.m_svo;
		//smo = pt.m_smo;									// smo from stress computation phase (from remodeled natural configuration)
		//sco = pt.m_sco;									// sco from stress computation phase (from remodeled natural configuration)
		//Fio = pt.m_Fio;
		//phic = pt.m_phic;								// converged phic from stress computation phase

		//phim = phimo/(J/Jo)*pow(J/Jo*phic/phico,eta);	// phim from <J*phim/phimo=(J*phic/phico)^eta>

		//F.right_polar(R, U); // mat3ds U; mat3d R;

		//U.eigen2(eigenval, eigenvec);
		////Ui = U.inverse();            					// Inverse of U
		////ui(Ui);

		// Compute current stresses (originally from Tangent function)
		sNm = smo;										// phim*smhato = phim*smo
		sNc = sco;										// phic*schato = phic*sco

		// Active
		//(Fio.inverse()).right_polar(Ro,Uo); //uo(Uo);	// Uo from polar decomposition
		//uo(Uo);

		//double lto = (Fio.inverse()*N[1]).norm();
		Sao = Tmax*(1.0-exp(-CB*CB))*(1.0-pow((lamM-1.0)/(lamM-lam0),2))*(lto*lto)*dyad(N[1]);
		sao = 1.0/Jo*(uo*(Sao*uo)).sym();

		//double lr = (F*(Fio*N[0])).norm();						// lr -> 1 for F -> Fo
		//double lt = (F*(Fio*N[1])).norm();						// lt -> 1 for F -> Fo

		//double rIrIo = ro/rIo*lt-(ro-rIo)/rIo*lr;				// rIrIo -> rIorIo = 1 for F -> Fo

		Cratio = CB-CS*(EPS*pow(rIrIo,-3)-1.0);
		sNa.zero(); // mat3ds sNa;
		if (Cratio>0) sNa = (1.0-exp(-Cratio*Cratio))/(1.0-exp(-CB*CB))*sao;

		// 2nd P-K stresses
		Se = (phieo*mu*Ge*Ge).sym();
		const mat3ds Sm = J*(ui*sNm*ui).sym();						// J*Ui*sNm*Ui
		const mat3ds Sc = J*(ui*sNc*ui).sym();						// J*Ui*sNc*Ui
		Sa = J*(ui*sNa*ui).sym();						// J*Ui*sNa*Ui

		Sf = phim*Sm+phic*Sc;							// J*Ui*sNf*Ui
		Sx = Se+Sf+phim*Sa;

		// Associated Cauchy stresses
		mat3ds sm = 1.0/J*(F*(Sm*F.transpose())).sym(); // const
		mat3ds sc = 1.0/J*(F*(Sc*F.transpose())).sym(); // const
		mat3ds sa = 1.0/J*(F*(Sa*F.transpose())).sym(); // const
		mat3ds sx = 1.0/J*(F*(Sx*F.transpose())).sym(); // const

		//p = 1.0/3.0/J*Sx.dotdot(C)-svo/(1.0-delta)*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0)-KfKi*inflam);	// Ups = 1 -> p

		// Compute tangent
		tens4dmm Ixsx = dyad1mm(I,sx); // const
		tens4dmm smxI = dyad1mm(sm,I); // const
		tens4dmm scxI = dyad1mm(sc,I); // const
		tens4dmm saxI = dyad1mm(sa,I); // const

		mat3ds tenr = dyad(F*(Fio*N[0])); // const						// Fio needed for consistency (from computation of lr)
		mat3ds tent = dyad(F*(Fio*N[1])); // const
		mat3ds tenz = dyad(F*(Fio*N[2])); // const

		tens4dmm Ixnrr = dyad1mm(I,tenr); // const
		tens4dmm Ixntt = dyad1mm(I,tent); // const

		// Contribution due to constant Cauchy stresses at constituent level
		tens4dmm cfss(0.0);

		sfpro.zero();
		sfpro(0,0) = eigenvec[0]*((phim*sNm+phic*sNc+phim*sNa)*eigenvec[0]);
		sfpro(1,1) = eigenvec[1]*((phim*sNm+phic*sNc+phim*sNa)*eigenvec[1]);
		sfpro(2,2) = eigenvec[2]*((phim*sNm+phic*sNc+phim*sNa)*eigenvec[2]);
		sfpro(0,1) = eigenvec[0]*((phim*sNm+phic*sNc+phim*sNa)*eigenvec[1]);
		sfpro(1,2) = eigenvec[1]*((phim*sNm+phic*sNc+phim*sNa)*eigenvec[2]);
		sfpro(0,2) = eigenvec[0]*((phim*sNm+phic*sNc+phim*sNa)*eigenvec[2]);

		vec3d Fxeigenvec[3];

		Fxeigenvec[0] = F*eigenvec[0];
		Fxeigenvec[1] = F*eigenvec[1];
		Fxeigenvec[2] = F*eigenvec[2];

		for (int i=0; i<3; i++) {

			mat3ds ten1 = dyad(Fxeigenvec[i]);

			for (int j=0; j<3; j++) {

				double component = sfpro(i,j) / pow(eigenval[i],3) / eigenval[j];

				mat3ds ten2 = dyads(Fxeigenvec[i],Fxeigenvec[j]);

				cfss -= component*dyad1mm(ten2,ten1);

				for (int k=0; k<3; k++) {

					if (k == i) continue;

					mat3ds ten3 = dyads(Fxeigenvec[j],Fxeigenvec[k]);
					mat3ds ten4 = dyads(Fxeigenvec[k],Fxeigenvec[i]);

					component = sfpro(i,j) / eigenval[i] / eigenval[j] / eigenval[k] / (eigenval[i] + eigenval[k]);

					cfss -= component*dyad1mm(ten3,ten4);
				}
			}
		}

		double dphiRm = phimo*eta*pow(J/Jo*phic/phico,eta-1.0)/(phimo*eta*pow(J/Jo*phic/phico,eta-1.0)+phico); // const
		double dphiRc = phico/(phimo*eta*pow(J/Jo*phic/phico,eta-1.0)+phico); // const

		cfss += dphiRm*(smxI+saxI) + dphiRc*scxI;

		// Contribution due to the ratio of vasocontrictors to vasodilators in the active stress
		tens4dmm saoxnrr = dyad1mm((R*sao*R.transpose()).sym(),tenr); // const
		tens4dmm saoxntt = dyad1mm((R*sao*R.transpose()).sym(),tent); // const

		// 1/J * FoF : [ J * phim * 1/(1.0-exp(-CB*CB)) * (Ui*sao*Ui) x d(1-exp(-Cratio^2))/d(C/2) ] : (Ft)o(Ft)
		tens4dmm cass = phim * 6.0*Cratio*CS*EPS*pow(rIrIo,-4)*exp(-Cratio*Cratio)/(1.0-exp(-CB*CB)) * (ro/rIo/lt*saoxntt-(ro-rIo)/rIo/lr*saoxnrr); // const

		// Contribution due to change in Cauchy stresses at constituent level (orientation only, for now)
		tens4dmm cpnss(0.0);

		//// Redundant code from Tangent
		//double lpo = (Fio.inverse()*Np).norm();					// original referential stretch for deposition stretch calculation
		//double lno = (Fio.inverse()*Nn).norm();					// idem for symmetric
		//
		//double lcp2 = (Gc*lpo)*(Gc*lpo);						// deposition stretch calculation (computational purposes)
		//double lcn2 = (Gc*lno)*(Gc*lno);						// idem for symmetric

		//double lz = (F*(Fio*N[2])).norm();						// lz -> 1 for F -> Fo

		//alpha = atan(tan(alpha)*pow(lt/lz,aexp));				// update alpha
		//Np = N[1]*sin(alpha)+N[2]*cos(alpha);					// update diagonal fiber vector
		//Nn = N[1]*sin(alpha)-N[2]*cos(alpha);					// idem for symmetric

		const double scphato = cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc);	// Constant stress magnitude at constituent level
		const double scnhato = cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc);

		const vec3d dNpdta = (N[1]-N[2]*tan(alpha))*pow(1+pow(tan(alpha),2),-1.5);	// d(Np)/d(tan(alpha))
		const vec3d dNndta = (N[1]+N[2]*tan(alpha))*pow(1+pow(tan(alpha),2),-1.5);

		const mat3ds ten1 = 1.0/Jo*dyads(R*(Uo*dNpdta),R*(Uo*Np));					// FoF : (Ui)o(Ui) : d(NpxNp)/d(tan(alpha)), with Jo and Uo needed for consistency (from computation of sco)
		const mat3ds ten2 = 1.0/Jo*dyads(R*(Uo*dNndta),R*(Uo*Nn));
		
		const mat3ds ten3 = aexp*tan(alpha)*(1.0/(lt*lt)*tent-1.0/(lz*lz)*tenz);		// 2*d(tan(alpha))/d(C) : (Ft)o(Ft)

		cpnss += (phic*betad) * scphato * dyad1mm(ten1,ten3);					// 1/J * FoF : [ J * phicp * scphato * (Ui)o(Ui) : 2*d(NpxNp)/d(C) ] : (Ft)o(Ft)
		cpnss += (phic*betad) * scnhato * dyad1mm(ten2,ten3);

		const tens4dmm cess = tens4dmm(ce);							// ce in tens4dmm form

		css = cess + cfss + cass + cpnss;

		//css += 1.0/3.0*(2.0*sx.tr()*IoIss-2.0*Ixsx-ddot(IxIss,css))
		//	 + svo/(1.0-delta)*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0)-KfKi*inflam)*(IxIss-2.0*IoIss)
		//	 - 3.0*svo/(1.0-delta)*KsKi*EPS*pow(rIrIo,-4)*(ro/rIo/lt*Ixntt-(ro-rIo)/rIo/lr*Ixnrr);
		css += 1.0/3.0*(2.0*sx.tr()*IoIss-2.0*Ixsx-ddot(IxIss,css))
			 + svo*(1.0 + adho*(sgr-1.0)/(partialtime-1.0))/(1.0-delta)*(1.0+KsKi*(EPS*pow(rIrIo,-3)-1.0)-KfKi*inflam)*(IxIss-2.0*IoIss)
			 - 3.0*svo*(1.0 + adho*(sgr-1.0)/(partialtime-1.0))/(1.0-delta)*KsKi*EPS*pow(rIrIo,-4)*(ro/rIo/lt*Ixntt-(ro-rIo)/rIo/lr*Ixnrr);

	}
	else {

		double    Jo = pt.m_Jo;
		double    Jh = pt.m_Jh;
		mat3d    Fio = pt.m_Fio;
		mat3d    Fih = pt.m_Fih;
		double phich = pt.m_phic;

		double phimh = phimo/(Jh/Jo)*pow(Jh/Jo*phich/phico,eta);	// Evolved homeostatic phimh from <Jh*phimh/phimo=(Jh*phich/phico)^eta>
		double phieh = phieo/(Jh/Jo);

		double Jdep = 0.9999;
		//double lm = 1.0e3*mu;
		double lm = bulkLM*mu;

		double lrh = (Fih.inverse()*(Fio*N[0])).norm();				// lrh -> 1 for Fh -> Fo
		double lth = (Fih.inverse()*(Fio*N[1])).norm();				// lth -> 1 for Fh -> Fo
		double lzh = (Fih.inverse()*(Fio*N[2])).norm();				// lzh -> 1 for Fh -> Fo

		alpha = atan(tan(alpha)*pow(lth/lzh,aexp));					// Remodeled (evolved homeostatic) alpha at h (alphah)
		Np = N[1]*sin(alpha)+N[2]*cos(alpha);						// Update diagonal fiber vector
		Nn = N[1]*sin(alpha)-N[2]*cos(alpha);						// idem for symmetric

		mat3ds Uo; mat3d Ro; (Fio.inverse()).right_polar(Ro,Uo);	// Uo from polar decomposition
		mat3ds Uh; mat3d Rh; (Fih.inverse()).right_polar(Rh,Uh);	// Uh from polar decomposition
		
		double lt = (F*(Uh.inverse()*(Uo*N[1]))).norm();
		double lz = (F*(Uh.inverse()*(Uo*N[2]))).norm();
		double lp = (F*(Uh.inverse()*(Uo*Np))).norm();
		double ln = (F*(Uh.inverse()*(Uo*Nn))).norm();
		
		double lmt2 = (Gm*lt)*(Gm*lt);
		double lct2 = (Gc*lt)*(Gc*lt);
		double lcz2 = (Gc*lz)*(Gc*lz);
		double lcp2 = (Gc*lp)*(Gc*lp);
		double lcn2 = (Gc*ln)*(Gc*ln);

		mat3ds tent = dyad(F*(Uh.inverse()*(Uo*N[1])));
		mat3ds tenz = dyad(F*(Uh.inverse()*(Uo*N[2])));
		mat3ds tenp = dyad(F*(Uh.inverse()*(Uo*Np)));
		mat3ds tenn = dyad(F*(Uh.inverse()*(Uo*Nn)));

		// Strain energy density (per unit current volume)
		double sede = phieh * ( mu/2.0*((F*Ge).dotdot(F*Ge)-3.0) );
		double sedm = phimh * ( cm/(4.0*dm)*(exp(dm*(lmt2-1.0)*(lmt2-1.0))-1.0));
		double sedc = phich * ( cc/(4.0*dc)*(exp(dc*(lct2-1.0)*(lct2-1.0))-1.0)*betat +
								cc/(4.0*dc)*(exp(dc*(lcz2-1.0)*(lcz2-1.0))-1.0)*betaz +
								cc/(4.0*dc)*(exp(dc*(lcp2-1.0)*(lcp2-1.0))-1.0)*betad +
								cc/(4.0*dc)*(exp(dc*(lcn2-1.0)*(lcn2-1.0))-1.0)*betad );
		
		et.m_v.y = phimh;
		et.m_v.z = phich;

		et.m_a.x = sede + sedm + sedc;
		et.m_a.y = phieo*(2.0*mu*(((F*Ge).transpose())*(F*Ge)).dotdot(dyad(N[1])))+
				   phimo*(2.0*cm*exp(dm*(lmt2-1.0)*(lmt2-1.0))*((lmt2-1.0)+(1.0+2.0*dm*(lmt2-1.0)*(lmt2-1.0))*pow(Gm,2))*pow(Gm,2))+
				   phico*(2.0*cc*exp(dc*(lct2-1.0)*(lct2-1.0))*((lct2-1.0)+(1.0+2.0*dc*(lct2-1.0)*(lct2-1.0))*pow(Gc,2))*pow(Gc,2)*betat+
						  2.0*cc*exp(dc*(lcp2-1.0)*(lcp2-1.0))*((lcp2-1.0)+(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,2)*pow(sin(alpha),2))*pow(Gc,2)*pow(sin(alpha),2)*betad+
						  2.0*cc*exp(dc*(lcn2-1.0)*(lcn2-1.0))*((lcn2-1.0)+(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,2)*pow(sin(alpha),2))*pow(Gc,2)*pow(sin(alpha),2)*betad);    // c_tttt
		et.m_a.z = phieo*(2.0*mu*(((F*Ge).transpose())*(F*Ge)).dotdot(dyad(N[2])))+
				   phico*(2.0*cc*exp(dc*(lcz2-1.0)*(lcz2-1.0))*((lcz2-1.0)+(1.0+2.0*dc*(lcz2-1.0)*(lcz2-1.0))*pow(Gc,2))*pow(Gc,2)*betaz+
						  2.0*cc*exp(dc*(lcp2-1.0)*(lcp2-1.0))*((lcp2-1.0)+(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,2)*pow(cos(alpha),2))*pow(Gc,2)*pow(cos(alpha),2)*betad+
						  2.0*cc*exp(dc*(lcn2-1.0)*(lcn2-1.0))*((lcn2-1.0)+(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,2)*pow(cos(alpha),2))*pow(Gc,2)*pow(cos(alpha),2)*betad);    // c_zzzz

		// Passive | consistent with: Sth = Jh*(Ui*sNt*Ui), sNt = phih/phio*sto, sto = 1.0/Jo*(Uo*(Sto*Uo)), Sto = phio*c*(lt2-1.0)*exp(d*(lt2-1.0)*(lt2-1.0))*(Gt*Gt)*dyad(N[1])
		mat3ds Se = (phieo*mu*Ge*Ge).sym();
		mat3ds Sm = (Jh/Jo*phimh)*(cm*(lmt2-1.0)*exp(dm*(lmt2-1.0)*(lmt2-1.0))*(Gm*Gm)*dyad(Uh.inverse()*(Uo*N[1])));
		mat3ds Sc = (Jh/Jo*phich)*(cc*(lct2-1.0)*exp(dc*(lct2-1.0)*(lct2-1.0))*(Gc*Gc)*dyad(Uh.inverse()*(Uo*N[1]))*betat +
								   cc*(lcz2-1.0)*exp(dc*(lcz2-1.0)*(lcz2-1.0))*(Gc*Gc)*dyad(Uh.inverse()*(Uo*N[2]))*betaz +
								   cc*(lcp2-1.0)*exp(dc*(lcp2-1.0)*(lcp2-1.0))*(Gc*Gc)*dyad(Uh.inverse()*(Uo* Np ))*betad +
								   cc*(lcn2-1.0)*exp(dc*(lcn2-1.0)*(lcn2-1.0))*(Gc*Gc)*dyad(Uh.inverse()*(Uo* Nn ))*betad );
		
		mat3ds Sx = Se + Sm + Sc;
		mat3ds sx = 1.0/J*((F*(Sx*F.transpose()))).sym();

		S = Sx + Ci*lm*log(Jdep*J/Jh);

		// Passive material stiffness
		tens4ds ce(0.0);
		tens4ds cf = (Jh/Jo*phimh)*(2.0*cm*(1.0+2.0*dm*(lmt2-1.0)*(lmt2-1.0))*exp(dm*(lmt2-1.0)*(lmt2-1.0))*pow(Gm,4)*dyad1s(tent))      +
					 (Jh/Jo*phich)*(2.0*cc*(1.0+2.0*dc*(lct2-1.0)*(lct2-1.0))*exp(dc*(lct2-1.0)*(lct2-1.0))*pow(Gc,4)*dyad1s(tent)*betat +
									2.0*cc*(1.0+2.0*dc*(lcz2-1.0)*(lcz2-1.0))*exp(dc*(lcz2-1.0)*(lcz2-1.0))*pow(Gc,4)*dyad1s(tenz)*betaz +
									2.0*cc*(1.0+2.0*dc*(lcp2-1.0)*(lcp2-1.0))*exp(dc*(lcp2-1.0)*(lcp2-1.0))*pow(Gc,4)*dyad1s(tenp)*betad +
									2.0*cc*(1.0+2.0*dc*(lcn2-1.0)*(lcn2-1.0))*exp(dc*(lcn2-1.0)*(lcn2-1.0))*pow(Gc,4)*dyad1s(tenn)*betad);
		cf /= J;

		tens4ds c = ce + cf;

		et.m_a.y = 2.0*sx.dotdot(dyad(F*(Uo*N[1])))/(F*(Uo*N[1])).norm2()+(c.dot(dyad(F*(Uo*N[1])))).dotdot(dyad(F*(Uo*N[1])))/pow((F*(Uo*N[1])).norm2(),2.0);
		et.m_a.z = 2.0*sx.dotdot(dyad(F*(Uo*N[2])))/(F*(Uo*N[2])).norm2()+(c.dot(dyad(F*(Uo*N[2])))).dotdot(dyad(F*(Uo*N[2])))/pow((F*(Uo*N[2])).norm2(),2.0);

		c += lm/J*(IxI-2.0*log(Jdep*J/Jh)*IoI);
		css = tens4dmm(c);
	}

	mat3ds s = 1.0/J*((F*(S*F.transpose()))).sym();
	stress = s;

	//et.m_v.x = sqrt(C.dotdot(dyad(N[0])));					// Radial stretch
	//et.m_v.y = sqrt(C.dotdot(dyad(N[1])));					// Circ. stretch
	//et.m_v.z = sqrt(C.dotdot(dyad(N[2])));					// Axial stretch
	//et.m_v.x = s.dotdot(dyad(F*N[0]))/(F*N[0]).norm2();			// Radial stress, just for plotting
	//et.m_v.y = s.dotdot(dyad(F*N[1]))/(F*N[1]).norm2();			// Circ.  stress, just for plotting
	//et.m_v.z = s.dotdot(dyad(F*N[2]))/(F*N[2]).norm2();			// Axial  stress, just for plotting
	double s_rr = s.dotdot(dyad(F*N[0]))/(F*N[0]).norm2();			// Radial stress, just for plotting
	double s_tt = s.dotdot(dyad(F*N[1]))/(F*N[1]).norm2();			// Circ.  stress, just for plotting
	double s_zz = s.dotdot(dyad(F*N[2]))/(F*N[2]).norm2();			// Axial  stress, just for plotting
	et.m_v.x = m_insult(mp)*(sgr-1.0)/(partialtime-1.0);        // Normalized insult value
	//et.m_v.y = inflam*KfKi;										// Inflammation
	//et.m_v.z = alpha/M_PI*180.;									// Diag. collagen angle (deg)
	//et.m_v.y = eta;
	//et.m_v.z = KsKi;
	//et.m_v.y = (F).dotdot(F);
	//et.m_v.z = mu;

	tangent = css;

	//pt.m_insult = m_delta;
	//if (abs(t-partialtime) <= eps) {
	if (t >= partialtime && trunc(t) == t) {
		//vec3d Xcur = et.m_rt;
		//feLog("%.6f %.6f %.6f 0.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", X.x, X.y, X.z, Xcur.x, Xcur.y, Xcur.z, et.m_a.x, et.m_a.y, et.m_a.z, et.m_v.x, et.m_v.y, et.m_v.z, et.m_J);
		//feLog("%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", X.x, X.y, X.z, et.m_a.x, et.m_a.y, et.m_a.z, et.m_v.x, et.m_v.y, et.m_v.z);
		feLog("%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", X.x, X.y, X.z, et.m_a.x, et.m_a.y, et.m_a.z, et.m_v.x, s_tt, s_zz, et.m_v.y, et.m_v.z);
	}
}
