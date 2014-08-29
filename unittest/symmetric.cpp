#include "pinocchio/spatial/fwd.hpp"
#include "pinocchio/spatial/se3.hpp"
#include "pinocchio/spatial/inertia.hpp"
#include "pinocchio/tools/timer.hpp"

#include <boost/random.hpp>
#include <assert.h>

#include "pinocchio/spatial/symmetric3.hpp"

#include <Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION( Eigen::Matrix3d );
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(se3::Symmetric3);

void timeSym3(const se3::Symmetric3 & S,
	      const se3::Symmetric3::Matrix3 & R,
	      se3::Symmetric3 & res)
{
  res = S.rotate(R);
}

void testSym3()
{
  using namespace se3;
  typedef Symmetric3::Matrix3 Matrix3;
  typedef Symmetric3::Vector3 Vector3;
  
  { 
    // op(Matrix3)
    {
      Matrix3 M = Matrix3::Random(); M = M*M.transpose();
      Symmetric3 S(M);
      assert( S.matrix().isApprox(M) );
    }

    // S += S
    {
      Symmetric3
	S = Symmetric3::Random(),
	S2 = Symmetric3::Random();
      Symmetric3 Scopy = S;
      S+=S2;
      assert( S.matrix().isApprox( S2.matrix()+Scopy.matrix()) );
    }

    // S + M
    {
      Symmetric3 S = Symmetric3::Random();
      Matrix3 M = Matrix3::Random(); M = M*M.transpose();

      Symmetric3 S2 = S + M;
      assert( S2.matrix().isApprox( S.matrix()+M ));

      S2 = S - M;
      assert( S2.matrix().isApprox( S.matrix()-M ));
    }

    // S*v
    {
      Symmetric3 S = Symmetric3::Random();
      Vector3 v = Vector3::Random(); 
      Vector3 Sv = S*v;
      assert( Sv.isApprox( S.matrix()*v ));
    }

    // Random
    for(int i=0;i<100;++i )
      {
	Matrix3 M = Matrix3::Random(); M = M*M.transpose();
	Symmetric3 S = Symmetric3::RandomPositive();
	Vector3 v = Vector3::Random();
	assert( (v.transpose()*(S*v))[0] > 0);
      }

    // Identity
    { 
      assert( Symmetric3::Identity().matrix().isApprox( Matrix3::Identity() ) );
    }

    // Skew2
    {
      Vector3 v = Vector3::Random();
      Symmetric3 vxvx = Symmetric3::SkewSq(v);

      Vector3 p = Vector3::UnitX();
      assert( (vxvx*p).isApprox( v.cross(v.cross(p)) ));
      p = Vector3::UnitY();
      assert( (vxvx*p).isApprox( v.cross(v.cross(p)) ));
      p = Vector3::UnitZ();
      assert( (vxvx*p).isApprox( v.cross(v.cross(p)) ));
    }

    // (i,j)
    {
	Matrix3 M = Matrix3::Random(); M = M*M.transpose();
	M << 1,2,4,2,3,5,4,5,6 ;
	Symmetric3 S(M); // = Symmetric3::RandomPositive();
	for(int i=0;i<3;++i)
	  for(int j=0;j<3;++j)
	    assert( S(i,j) == M(i,j) );
    }
  }

  // SRS
  {
    Symmetric3 S = Symmetric3::RandomPositive();
    Matrix3 R = (Eigen::Quaterniond(Eigen::Matrix<double,4,1>::Random())).normalized().matrix();
    
    Symmetric3 RSRt = S.rotate(R);
    assert( RSRt.matrix().isApprox( R*S.matrix()*R.transpose() ));

    Symmetric3 RtSR = S.rotate(R.transpose());
    assert( RtSR.matrix().isApprox( R.transpose()*S.matrix()*R ));

  }

  // Time test 
  {
    const int NBT = 100000;
    Symmetric3 S = Symmetric3::RandomPositive();

    std::vector<Symmetric3> Sres (NBT);
    std::vector<Matrix3> Rs (NBT);
    for(int i=0;i<NBT;++i) 
      Rs[i] = (Eigen::Quaterniond(Eigen::Matrix<double,4,1>::Random())).normalized().matrix();

    StackTicToc timer(StackTicToc::US); timer.tic();
    SMOOTH(NBT)
      {
	timeSym3(S,Rs[_smooth],Sres[_smooth]);
      }
    timer.toc(std::cout,NBT);
  }
}
/* --- METAPOD ---------------------------------------------- */
#include <metapod/tools/spatial/lti.hh>
#include <metapod/tools/spatial/rm-general.hh>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(metapod::Spatial::ltI<double>);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(metapod::Spatial::RotationMatrixTpl<double>);

void timeLTI(const metapod::Spatial::ltI<double>& S,
	     const metapod::Spatial::RotationMatrixTpl<double>& R, 
	     metapod::Spatial::ltI<double> & res)
{
  res = R.rotSymmetricMatrix(S);
}

void testLTI()
{
  using namespace metapod::Spatial;

  typedef ltI<double> Sym3;
  typedef Eigen::Matrix3d Matrix3;
  typedef RotationMatrixTpl<double> R3;

  Matrix3 M = Matrix3::Random();
  Sym3 S(M),S2;

  R3 R; R.randomInit();

  R.rotTSymmetricMatrix(S);
  timeLTI(S,R,S2);

  assert( S2.toMatrix().isApprox( R.toMatrix().transpose()*S.toMatrix()*R.toMatrix()) );
  
  const int NBT = 100000;
  std::vector<Sym3> Sres (NBT);
  std::vector<R3> Rs (NBT);
  for(int i=0;i<NBT;++i) 
    Rs[i].randomInit();
  
  StackTicToc timer(StackTicToc::US); timer.tic();
  SMOOTH(NBT)
    {
      timeLTI(S,Rs[_smooth],Sres[_smooth]);
    }
  timer.toc(std::cout,NBT);
  //std::cout << Rs[std::rand() % NBT] << std::endl;
  
}

void timeSelfAdj( const Eigen::Matrix3d & A,
		  const Eigen::Matrix3d & Sdense,
		  Eigen::Matrix3d & ASA )
{
  typedef se3::Inertia::Symmetric3 Sym3;
  Sym3 S(Sdense);
  ASA.triangularView<Eigen::Upper>()
    = A * S * A.transpose();
}

void testSelfAdj()
{
  using namespace se3;
  typedef Inertia::Matrix3 Matrix3;
  typedef Inertia::Symmetric3 Sym3;

  Matrix3 M = Inertia::Matrix3::Random();
  Sym3 S(M);

  Matrix3 M2 = Inertia::Matrix3::Random();
  M.triangularView<Eigen::Upper>() = M2;

  Matrix3 A = Matrix3::Random(), ASA1, ASA2;
  ASA1.triangularView<Eigen::Upper>() = A * S * A.transpose();
  timeSelfAdj(A,M,ASA2);
  assert(ASA1.isApprox(ASA2));

  StackTicToc timer(StackTicToc::US); timer.tic();
  SMOOTH(1000000)
    {
      timeSelfAdj(A,M,ASA2);
    }
  timer.toc(std::cout,1000000);
}


int main()
{
  //testSelfAdj();
  testLTI();
  testSym3();
  testSym3();
  testLTI();
  
  std::cout << std::endl;
  return 0;
}