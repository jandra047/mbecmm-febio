#pragma once
//=============================================================================
// This plugin illustrates how to create a mechanobiologically equilibrated constrained mixture model. 
// It requires FEBio 3.2 (or up)
//
// Authors: Martin Pfaller, Marcos Latorre, David Li
// Copyright (c) 2021
// All rights reserved
//
//=============================================================================

//-----------------------------------------------------------------------------
// We need to include this file since our new material class will inherit from
// FEElasticMaterial, which is defined in this include file.
#include "FEBioMech/FEElasticMaterial.h"

class GRMaterialPoint : public FEMaterialPoint
{
public:
	GRMaterialPoint(FEMaterialPoint *pt);

	FEMaterialPoint* Copy() override;

	void Init() override;

	void Serialize(DumpStream& ar) override;

public:
	// Original (o) homeostatic data
	double		m_Jo;		//!< Jacobian at o
	double		m_svo;		//!< Volumetric stress at o
	mat3ds		m_smo;		//!< Cauchy stress tensor for smooth muscle cells at o
	mat3ds		m_sco;		//!< Cauchy stress tensor for all collagen fiber families at o
	mat3d		m_Fio;		//!< Inverse of deformation gradient tensor at o
	double		m_Jh;		//!< Jacobian at h
	mat3d		m_Fih;		//!< Inverse of deformation gradient tensor at h

	// Evolved homeostatic (h) data
	double		m_phic;		//!< Total mass fraction of all collagen fiber families at h
	double		m_Iemax;	//!< Maximum value of Ie achieved over the loading history up until the current G&R time
};

//-----------------------------------------------------------------------------
// This material class implements a mechanobiologically equilibrated constrained mixture model. 
// Since it contains a hyperelatic material step, it needs to inherit from FEElasticMaterial.
class FEMbeCmm : public FEElasticMaterial
{
public:
	// Returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() override { return new GRMaterialPoint(new FEElasticMaterialPoint); }

public:
	// The constructor is called when an instance of this class is created.
	// All classes registered by the framework must take the FEModel* as the only
	// parameter in the constructor, even if the class does not need it (which most often
	// will be the case). For material classes, the FEModel parameter is passed to the 
	// base class in the initialization list.

	// Setting m_secant_tangent = true so FESolidMaterial uses SecantTangent
	// (allows minor symmetry only tangents) instead of Tangent (minor and major symmetries)
	FEMbeCmm(FEModel* pfem) : FEElasticMaterial(pfem) {
		m_secant_tangent = true;
	}

public:
	// Simulation
	double rIo;         // Initial inner radius
	double lo;          // Initial axial length
	double endtime;     // Simulation end time
	double partialtime; // G&R end time, usually = endtime
	double Jdep;        // Near incompressibility
	double bulkLM;      // Incompressibility multiplier
	double imper;       // Tortuosity amount
	double hwaves;      // Tortuosity period

	// Structure
	double phieo;      // Initial elastin mass fraction
	double phimo;      // Initial muscle mass fraction
	double phico;      // Initial collagen mass fraction
	double betat;      // Initial circumferential collagen fraction
	double betaz;      // Initial axial collagen fraction
	//double betad;      // Initial diagonal collagen fraction: 0.5*(1.0 - betat - betaz)
	double alphao;      // Initial diagonal collagen orientation (rad)

	// Material properties
	double mub;        // Elastin shear modulus, baseline
	double cm;         // Muscle modulus a
	double dm;         // Muscle modulus b
	double ccb;        // Collagen modulus a, baseline
	double dc;         // Collagen modulus b

	double Get;        // Elastin deposition stretch circumferential
	double Gez;        // Elastin deposition stretch axial
	double Gm;         // Muscle deposition stretch
	double Gcb;        // Collagen deposition stretch, baseline

	// Active
	double Tmaxb;      // Muscle tone, baseline
	double CB;         // Shape parameter s.t. (1-exp(-CB^2)) = 0.5
	//double CS;         // Shape parameter s.t. (1-exp(-C^2)) = 0.0 for lt = 1/(1+CB/CS)^(1/3) = 0.7 and (1-exp(-C^2)) = 0.75 for lt = 2.0
	double lamM;       // Maximal contraction stretch
	double lam0;       // Minimal contraction stretch

	// G&R
	double aexp;       // Allow for collagen reorientation
	double etab;       // Combined production-removal, baseline
	double KsKib;      // Shear-to-intramural gain ratio, baseline
	double infl;       // Inflammation coefficient
	double KfKicen;    // Inflam-to-intramural gain ratio
	double deltab;     // Mechanosensing, baseline
	//double ksi;        // Endothelial health (not implemented)

	// Aneurysm
	double asym;	   // Allow for asymmetric insult
	double zod;		   // Axial width (width = lo/2/zod)
	double nuz;		   // Axial slope
	double tod;		   // Circumferential width (width = pi/tod)
	double nut;		   // Circumferential slope
	double zo;         // Aneurysm apex axial location
	double to;         // Aneurysm apex circumferential location (*pi)
	
	// Insult
	double insmu;      // % reduction elastic fiber integrity
	double inscc;      // % reduction collagen cross-linking
	double insTmax;    // % reduction contractility
	double insdelta;   // Dysfunctional mechanosensing
	double insGc;      // Dysfunctional mechanoregulation
	double agimu;      // Uniform reduction elastic fiber integrity with aging
	double etacen;     // Combined production-removal in aneurysm center
	double KsKicen;    // Shear-to-intramural gain ratio in aneurysm center
	
	// Dev
	double adho;       // Percent change in homeostatic intramural stress (adaptive homeostasis)
	FEParamDouble m_insult; // Normalized insult value mapped from MeshData

	// Output Directory options
	//D:\Yale\FEBio\NumKnock\plugins\3.2.0\
	//..\..\..\NumKnock\plugins\3.2.0\
	//$(SolutionDir)$(Platform)\$(Configuration)\

public:
	void StressTangent(FEMaterialPoint& mp, mat3ds& stress, tens4dmm& tangent);

	// This function calculates the spatial (i.e., Cauchy or true) stress.
	// It takes one parameter, the FEMaterialPoint, and returns a mat3ds object,
	// which is a symmetric second-order tensor.
	virtual mat3ds Stress(FEMaterialPoint& pt){
		mat3ds stress;
		tens4dmm tangent;
		StressTangent(pt, stress, tangent);
		return stress;
	}

	// This function calculates the spatial elasticity tangent tensor. 
	// It takes one parameter, the FEMaterialPoint, and returns a tens4ds object,
	// which is a fourth-order tensor with major and minor symmetries.

	// Added by David Li (Sept 2021) to circumvent Error C4716 FEMbeCMM::Tangent must return a value
	//virtual tens4ds Tangent(FEMaterialPoint& pt) {};
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	// Minor symmetries only
	virtual tens4dmm SecantTangent(FEMaterialPoint& pt) {
		mat3ds stress;
		tens4dmm tangent;
		StressTangent(pt, stress, tangent);
		return tangent;
	}

	// This macro defines that the class will define a material parameter list.
	// The material parameter list itself is defined elsewhere (i.e., in the .cpp file.)
	DECLARE_FECORE_CLASS();
};
