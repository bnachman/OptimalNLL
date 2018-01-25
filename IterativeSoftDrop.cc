#include <sstream>
#include "IterativeSoftDrop.h"


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{



  // *************************************************************************************
  //
  // Execute the Iterative Soft Drop Algorithm
  //
  // *************************************************************************************



  // internal function called recursively by the CA variant of getSubjets (see below)
  void IterativeSoftDrop::_FindHardPts(const fastjet::PseudoJet & this_jet, 
				       std::vector<double> & t_parts) const
  {
    fastjet::PseudoJet parent1(0.0,0.0,0.0,0.0), parent2(0.0,0.0,0.0,0.0);
    bool had_parents=this_jet.validated_cs()->has_parents(this_jet,parent1,parent2);

    if ( !had_parents )
      {
	t_parts.push_back(-1.0);  return;
      } // stop recursion on this branch

    if (had_parents && parent1.plain_distance(parent2) < (_theta_cut*_theta_cut) )
      {
	t_parts.push_back(-2.0);  return;
      } // stop recursion on this branch

    if (parent1.perp() < parent2.perp()) std::swap(parent1,parent2);

    double pt1=parent1.perp();
    double pt2=parent2.perp();
    double totalpt=pt1+pt2;

    double zij = pt2/totalpt;
    double softdropcondition = _z_cut * pow( (sqrt(parent1.plain_distance(parent2)) / _Ro) , _beta );

    if ( zij > softdropcondition)   
      {
	_FindHardPts(parent1, t_parts);
	t_parts.push_back(zij);
      } // continue recursion on hard branch and do add fraction

    else{
      _FindHardPts(parent1, t_parts);
      // continue recursion on hard branch and don't add fraction
    }
    return;
  }

  void IterativeSoftDrop::_FindHardAngles(const fastjet::PseudoJet & this_jet, 
					  std::vector<double> & t_parts) const
  {
    fastjet::PseudoJet parent1(0.0,0.0,0.0,0.0), parent2(0.0,0.0,0.0,0.0);
    bool had_parents=this_jet.validated_cs()->has_parents(this_jet,parent1,parent2);

    if ( !had_parents )
      {
	t_parts.push_back(-1.0);  return;
      } // stop recursion on this branch

    if (had_parents && parent1.plain_distance(parent2) < (_theta_cut*_theta_cut) )
      {
	t_parts.push_back(-2.0);  return;
      } // stop recursion on this branch

    if (parent1.perp() < parent2.perp()) std::swap(parent1,parent2);

    double pt1=parent1.perp();
    double pt2=parent2.perp();
    double totalpt=pt1+pt2;

    double zij = pt2/totalpt;
    double softdropcondition = _z_cut * pow( (sqrt(parent1.plain_distance(parent2)) / _Ro) , _beta );

    if ( zij > softdropcondition)   
      {
	_FindHardAngles(parent1, t_parts);
	t_parts.push_back( sqrt(parent1.plain_distance(parent2)) );
      } // continue recursion on hard branch and do add angle

    else{
      _FindHardAngles(parent1, t_parts);
      // continue recursion on hard branch and don't add angle
    }
    return;
  }


  // the actual n_CA function, which is used internally; 
  // most of the hard work done by _FindHardPts which is called by getPtFracss 
  unsigned int IterativeSoftDrop::_n_CA(const fastjet::PseudoJet &jet) const
  {
    return (unsigned int)(IterativeSoftDrop::getPtFracs(jet).size() - 1);
  }

  // get the Pt fractions
  std::vector<double> IterativeSoftDrop::getPtFracs(const PseudoJet& jet) const
  {
    //JetAlgorithm algorithm = cambridge_algorithm;
    //double jet_rad = fastjet::JetDefinition::max_allowable_R; 
    // want to make sure we capture all the constituents in a single jet
    //JetDefinition jetDef = JetDefinition(algorithm, jet_rad, E_scheme, Best);
    //ClusterSequence clust_seq(jet.constituents(), jetDef);
    //std::vector<PseudoJet> ca_jets = sorted_by_pt(clust_seq.inclusive_jets());
 
    std::vector<double> fracs;
    _FindHardPts(jet, fracs);


    return fracs;
  }

  // get the Branching Angles
  std::vector<double> IterativeSoftDrop::getBranchAngles(const PseudoJet& jet) const
  {
    //JetAlgorithm algorithm = cambridge_algorithm;
    //double jet_rad = fastjet::JetDefinition::max_allowable_R; 
    // want to make sure we capture all the constituents in a single jet
    //JetDefinition jetDef = JetDefinition(algorithm, jet_rad, E_scheme, Best);
    //ClusterSequence clust_seq(jet.constituents(), jetDef);
    //std::vector<PseudoJet> ca_jets = sorted_by_pt(clust_seq.inclusive_jets());
 
    std::vector<double> angles;
    _FindHardAngles(jet, angles);


    return angles;
  }

  //-------------------------------
  unsigned int IterativeSoftDrop::result(const PseudoJet& jet) const {
   
    // if jet does not have constituents, throw error
    if (!jet.has_constituents()) throw Error("IterativeSoftDrop called on jet with no constituents.");
   
    return _n_CA(jet);
  }


  //---------------------------------
  std::string IterativeSoftDrop::description() const {
    std::ostringstream oss;
    oss << "IterativeSoftDrop using ";
    oss << "parameters z_cut = " << _z_cut << 
      ", beta = " << _beta << ", theta_cut = " << _theta_cut; 
    return oss.str();
  }


} // namespace contrib

FASTJET_END_NAMESPACE
