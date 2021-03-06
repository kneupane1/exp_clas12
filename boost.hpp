#ifndef BOOST_H_GUARD
#define BOOST_H_GUARD

#include "TLorentzVector.h"
#include "constants.hpp"
#include "physics.hpp"

class Boost_cm {
private:
  TLorentzVector *_beam;
  TLorentzVector *_elec;
  TLorentzVector *_target;
  TLorentzVector *_prot;
  TLorentzVector *_pip;
  TLorentzVector *_pim;

  TLorentzVector *_q_cm;
  TLorentzVector *_p_mu_prime_cm;
  TLorentzVector *_pip_mu_prime_cm;
  TLorentzVector *_pim_mu_prime_cm;

  bool _hasE;
  bool _hasP;
  bool _hasPip;
  bool _hasPim;

public:
  Reaction();
  Reaction(TLorentzVector *beam);
  ~Reaction();

  TLorentzVector e_mu_prime(); // maile thapeko
  TLorentzVector p_mu_prime();
  TLorentzVector pip_mu_prime();
  TLorentzVector pim_mu_prime();
  //  TLorentzVector kp_mu_prime();
  // TLorentzVector km_mu_prime();
  TLorentzVector q_cm(); // maile thapeko
  TLorentzVector p_mu_prime_cm();
  TLorentzVector pip_mu_prime_cm();
  TLorentzVector pim_mu_prime_cm();
  double theta_();
  void PhiAng(float phi);
  void AlphaCalc();

  //  void boost_fn(/*TLorentzVector four_vect, TLorentzVector e_mu,
  // TLorentzVector e_mu_prime);

  bool twoPionEvent();
};

#endif
