#ifndef UNAP_SMOOTHER_HPP
#define UNAP_SMOOTHER_HPP
namespace UNAP {
/// \brief abstract base-class for matrix Smoothers
class Smoother {
 protected:
  /// \brief communication domain
  Communicator *commcator_;

 public:
  /// \brief constructor
  /// \param other_comm communication domain
  Smoother(Communicator *other_comm) : commcator_(other_comm){};

  /// \brief destructor
  virtual ~Smoother() {}

  /// \brief smooth the unknown vector x
  /// \param x the unknown vector x
  /// \param A Matrix
  /// \param b rhs field
  /// \param nSweeps number of sweeps
  virtual void smooth(scalarVector &x, const Matrix &A, const scalarVector &b,
                      const label nSweeps) const = 0;

  /// \brief initialize the Smoother
  virtual void init() const = 0;

  /// \brief set communicator
  /// \param other_comm communication domain
  void setCommunicator(Communicator *other_comm) { commcator_ = other_comm; }

  /// \brief get communicator
  /// \param other_comm communication domain
  Communicator *getCommunicator() const { return commcator_; }
};
}  // namespace UNAP
#endif
