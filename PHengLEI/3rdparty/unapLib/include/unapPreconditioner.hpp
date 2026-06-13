#ifndef UNAP_PRECONDITIONER_HPP
#define UNAP_PRECONDITIONER_HPP

#include "unapVector.hpp"

namespace UNAP {
/// \brief abstract base-class for matrix Preconditioners
class Preconditioner {
 private:
  /// \brief brief description
  scalarVector *rDTempPtr_;

 protected:
  /// \brief communication domain
  Communicator *commcator_;

 public:
  /// \brief constructor
  /// \param other_comm communication domain
  Preconditioner(Communicator *other_comm) : commcator_(other_comm){};

  /// \brief destructor
  virtual ~Preconditioner() {}

  /// \brief return wA the preconditioned form of residual rA, \f$ w=P^{-1} r
  /// \f$, and \f$ P^{-1}Ax=P^{-1}b \f$ \param wA resulting vector,
  /// preconditioned by rA \param rA Preconditioner vector
  virtual void precondition(scalarVector &wA, const scalarVector &rA) const;

  /// \brief return the tempororary vector
  virtual const scalarVector &rD() const { return *rDTempPtr_; }

  /// \brief set communicator
  /// \param other_comm communication domain
  void setCommunicator(Communicator *other_comm) { commcator_ = other_comm; }

  /// \brief get communicator
  Communicator *getCommunicator() const { return commcator_; }
};
}  // namespace UNAP

#endif
