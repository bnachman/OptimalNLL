#ifndef __FASTJET_CONTRIB_ITERATIVESOFTDROP_H__
#define __FASTJET_CONTRIB_ITERATIVESOFTDROP_H__

#include <fastjet/internal/base.hh>
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include <vector>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{


  class IterativeSoftDrop : public FunctionOfPseudoJet<unsigned int> {
    // *************************************************************************************
    //
    // constructor Iterative SoftDrop
    //  https://arxiv.org/pdf/1704.06266.pdf
    //
    // ~ input:  z_cut is part of the soft drop criterion 
    // ~ input:  beta is exponent of deltaR/R 
    // ~ input:  Theta_cut is the branching angle cut
    //
    // ~ return: n_sd the integer number of branchings passing the cuts 
    // ~ return: z_n a vector of the transverse momentum fractions 
    // ~ return: Theta_n a vector of the branching angles
    //
    //   see original reference for details and discussion
    //
    // *************************************************************************************

  public:
  IterativeSoftDrop( double z_cut, double beta, double theta_cut, double Ro ) : 
    _z_cut(z_cut), _beta(beta), _theta_cut(theta_cut), _Ro(Ro) {}

    /// default dtor
    virtual ~IterativeSoftDrop(){}

    /// returns the number of branchings passing the cuts 
    unsigned int result(const PseudoJet& jet) const;

    /// get the Pt fraction vector
    std::vector<double> getPtFracs(const fastjet::PseudoJet &jet) const;
    /// get the branching angle vector
    std::vector<double> getBranchAngles(const fastjet::PseudoJet &jet) const;

    /// returns the description of the class
    std::string description() const; 

  private: 
    double _z_cut, _beta, _theta_cut, _Ro; // parameters defining n_CA

    // internal function that calculates the result
    unsigned int _n_CA(const fastjet::PseudoJet &jet) const;

    // internal function called recursively by _n_CA
    void _FindHardPts(const fastjet::PseudoJet & this_jet, 
		      std::vector<double> & t_parts) const;

    // internal function called recursively by _n_CA
    void _FindHardAngles(const fastjet::PseudoJet & this_jet, 
			 std::vector<double> & t_parts) const;


  }; // end class IterativeSoftDrop

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_ITERATIVESOFTDROP_H__
