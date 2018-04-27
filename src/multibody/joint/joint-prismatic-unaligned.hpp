//
// Copyright (c) 2015-2017 CNRS
// Copyright (c) 2016 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//
// This file is part of Pinocchio
// Pinocchio is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// Pinocchio is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// Pinocchio If not, see
// <http://www.gnu.org/licenses/>.

#ifndef __se3_joint_prismatic_unaligned_hpp__
#define __se3_joint_prismatic_unaligned_hpp__

#include "pinocchio/macros.hpp"
#include "pinocchio/multibody/joint/joint-base.hpp"
#include "pinocchio/multibody/constraint.hpp"
#include "pinocchio/spatial/inertia.hpp"

#include <stdexcept>

namespace se3
{

  struct MotionPrismaticUnaligned;
  
  template <>
  struct traits <MotionPrismaticUnaligned>
  {
    typedef double Scalar;
    typedef Eigen::Matrix<Scalar,3,1,0> Vector3;
    typedef Eigen::Matrix<Scalar,6,1,0> Vector6;
    typedef Eigen::Matrix<Scalar,6,6,0> Matrix6;
    typedef EIGEN_REF_CONSTTYPE(Vector6) ToVectorConstReturnType;
    typedef EIGEN_REF_TYPE(Vector6) ToVectorReturnType;
    typedef Vector3 AngularType;
    typedef Vector3 LinearType;
    typedef const Vector3 ConstAngularType;
    typedef const Vector3 ConstLinearType;
    typedef Matrix6 ActionMatrixType;
    typedef MotionTpl<Scalar,0> MotionPlain;
    enum {
      LINEAR = 0,
      ANGULAR = 3
    };
  }; // traits MotionPrismaticUnaligned

  struct MotionPrismaticUnaligned : MotionBase <MotionPrismaticUnaligned>
  {
    MOTION_TYPEDEF(MotionPrismaticUnaligned);

    MotionPrismaticUnaligned () : axis(Vector3::Constant(NAN)), rate(NAN) {}
    MotionPrismaticUnaligned (const Vector3 & axis, const Scalar rate) : axis(axis), rate(rate) {}

    Vector3 axis;
    Scalar rate;

    operator Motion() const { return Motion(axis*rate, Vector3::Zero());}
    
    template<typename Derived>
    void addTo(MotionDense<Derived> & v) const
    {
      v.linear() += axis * rate;
    }
  }; // struct MotionPrismaticUnaligned

  inline const MotionPrismaticUnaligned & operator+ (const MotionPrismaticUnaligned & m, const BiasZero &)
  { return m; }

  inline Motion operator+ (const MotionPrismaticUnaligned & m1, const Motion & m2)
  {
    return Motion(m1.rate*m1.axis + m2.linear(), m2.angular());
  }

  struct ConstraintPrismaticUnaligned;
  template <>
  struct traits <ConstraintPrismaticUnaligned>
  {
    typedef double Scalar;
    typedef Eigen::Matrix<Scalar,3,1,0> Vector3;
    typedef Eigen::Matrix<Scalar,4,1,0> Vector4;
    typedef Eigen::Matrix<Scalar,6,1,0> Vector6;
    typedef Eigen::Matrix<Scalar,3,3,0> Matrix3;
    typedef Eigen::Matrix<Scalar,4,4,0> Matrix4;
    typedef Eigen::Matrix<Scalar,6,6,0> Matrix6;
    typedef Matrix3 Angular_t;
    typedef Vector3 Linear_t;
    typedef const Matrix3 ConstAngular_t;
    typedef const Vector3 ConstLinear_t;
    typedef Matrix6 ActionMatrix_t;
    typedef Eigen::Quaternion<Scalar,0> Quaternion_t;
    typedef SE3Tpl<Scalar,0> SE3;
    typedef ForceTpl<Scalar,0> Force;
    typedef MotionTpl<Scalar,0> Motion;
    typedef Symmetric3Tpl<Scalar,0> Symmetric3;
    enum {
      LINEAR = 0,
      ANGULAR = 3
    };
    typedef Eigen::Matrix<Scalar,1,1,0> JointMotion;
    typedef Eigen::Matrix<Scalar,1,1,0> JointForce;
    typedef Eigen::Matrix<Scalar,6,1> DenseBase;
    typedef DenseBase MatrixReturnType;
    typedef const DenseBase ConstMatrixReturnType;
  }; // traits ConstraintPrismaticUnaligned

    struct ConstraintPrismaticUnaligned : ConstraintBase <ConstraintPrismaticUnaligned>
    {
      SPATIAL_TYPEDEF_NO_TEMPLATE(ConstraintPrismaticUnaligned);
      enum { NV = 1, Options = 0 };
      typedef traits<ConstraintPrismaticUnaligned>::JointMotion JointMotion;
      typedef traits<ConstraintPrismaticUnaligned>::JointForce JointForce;
      typedef traits<ConstraintPrismaticUnaligned>::DenseBase DenseBase;
      
      ConstraintPrismaticUnaligned () : axis (Vector3::Constant(NAN)) {}
      ConstraintPrismaticUnaligned (const Vector3 & axis) : axis(axis) {}
               
      Vector3 axis;

      template<typename D>
      MotionPrismaticUnaligned operator* (const Eigen::MatrixBase<D> & v) const
      { 
      	EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(D,1);
      	return MotionPrismaticUnaligned(axis,v[0]); 
      }

      Vector6 se3Action (const SE3 & m) const
      {
        Vector6 res;
        res.head<3> () = m.rotation()*axis;
        res.tail<3>().setZero();
        return res;
      }

      int nv_impl() const { return NV; }

      struct TransposeConst
      {
        typedef traits<ConstraintPrismaticUnaligned>::Scalar Scalar;
        typedef traits<ConstraintPrismaticUnaligned>::Force Force;
        typedef traits<ConstraintPrismaticUnaligned>::Vector6 Vector6;
        
      	const ConstraintPrismaticUnaligned & ref; 
      	TransposeConst(const ConstraintPrismaticUnaligned & ref) : ref(ref) {}
        
        template<typename Derived>
        Eigen::Matrix<
        typename EIGEN_DOT_PRODUCT_RETURN_TYPE(Vector3,typename ForceDense<Derived>::ConstLinearType),
        1,1>
      	operator* (const ForceDense<Derived> & f) const
      	{
          typedef Eigen::Matrix<
          typename EIGEN_DOT_PRODUCT_RETURN_TYPE(Vector3,typename ForceDense<Derived>::ConstLinearType),
          1,1> ReturnType;
          
          ReturnType res; res[0] = ref.axis.dot(f.linear());
          return res;
      	}

        /* [CRBA]  MatrixBase operator* (Constraint::Transpose S, ForceSet::Block) */
        template<typename D>
        friend
#if EIGEN_VERSION_AT_LEAST(3,2,90)
        const Eigen::Product<
        Eigen::Transpose<const Vector3>,
        typename Eigen::MatrixBase<const D>::template NRowsBlockXpr<3>::Type
        >
#else
        const typename Eigen::ProductReturnType<
        Eigen::Transpose<const Vector3>,
        typename Eigen::MatrixBase<const D>::template NRowsBlockXpr<3>::Type
        >::Type
#endif
        operator* (const TransposeConst & tc, const Eigen::MatrixBase<D> & F)
        {
          EIGEN_STATIC_ASSERT(D::RowsAtCompileTime==6,THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE)
          /* Return ax.T * F[1:3,:] */
          return tc.ref.axis.transpose () * F.template topRows<3> ();
        }

      };
      TransposeConst transpose() const { return TransposeConst(*this); }


      /* CRBA joint operators
       *   - ForceSet::Block = ForceSet
       *   - ForceSet operator* (Inertia Y,Constraint S)
       *   - MatrixBase operator* (Constraint::Transpose S, ForceSet::Block)
       *   - SE3::act(ForceSet::Block)
       */
      DenseBase matrix_impl() const
      {
        DenseBase S;
      	S << axis, Vector3::Zero();
      	return S;
      }
      
      DenseBase motionAction(const Motion & m) const
      {
        DenseBase res;
        res << m.angular().cross(axis), Vector3::Zero();
        
        return res;
      }
      
    }; // struct ConstraintPrismaticUnaligned


    inline Motion operator^ (const Motion & m1, const MotionPrismaticUnaligned & m2)
    {
      /* m1xm2 = [ v1xw2 + w1xv2; w1xw2 ] = [ v1xw2; w1xw2 ] */
      const Motion::Vector3 & w1 = m1.angular();
      const Motion::Vector3 & v2 = m2.axis * m2.rate;
      return Motion (w1.cross(v2), Motion::Vector3::Zero());
    }

    /* [CRBA] ForceSet operator* (Inertia Y,Constraint S) */
    inline Eigen::Matrix<double,6,1>
    operator* (const Inertia & Y, const ConstraintPrismaticUnaligned & cpu)
    { 
      /* YS = [ m -mcx ; mcx I-mcxcx ] [ v ; 0 ] = [ mv ; mcxv ] */
      const double & m                = Y.mass();
      const Inertia::Vector3 & c      = Y.lever();

      Eigen::Matrix<double,6,1> res;
      res.head<3>() = m*cpu.axis;
      res.tail<3>() = c.cross(res.head<3>());
      return res;
    }
  
  /* [ABA] Y*S operator (Inertia Y,Constraint S) */
  inline
#if EIGEN_VERSION_AT_LEAST(3,2,90)
  const Eigen::Product<
  Eigen::Block<const Inertia::Matrix6,6,3>,
  ConstraintPrismaticUnaligned::Vector3,
  Eigen::DefaultProduct>
#else
  const Eigen::ProductReturnType<
  Eigen::Block<const Inertia::Matrix6,6,3>,
  const ConstraintPrismaticUnaligned::Vector3
  >::Type
#endif
  operator*(const Inertia::Matrix6 & Y, const ConstraintPrismaticUnaligned & cpu)
  {
    return Y.block<6,3> (0,Inertia::LINEAR) * cpu.axis;
  }

  
    namespace internal 
    {
      template<>
      struct SE3GroupAction<ConstraintPrismaticUnaligned >  
      { typedef Eigen::Matrix<double,6,1> ReturnType; };
      
      template<>
      struct MotionAlgebraAction<ConstraintPrismaticUnaligned>
      { typedef Eigen::Matrix<double,6,1> ReturnType; };
    }

    struct JointPrismaticUnaligned;
    template<>
    struct traits<JointPrismaticUnaligned>
    {
      enum {
        NQ = 1,
        NV = 1
      };
      typedef double Scalar;
      typedef JointDataPrismaticUnaligned JointDataDerived;
      typedef JointModelPrismaticUnaligned JointModelDerived;
      typedef ConstraintPrismaticUnaligned Constraint_t;
      typedef SE3 Transformation_t;
      typedef MotionPrismaticUnaligned Motion_t;
      typedef BiasZero Bias_t;
      typedef Eigen::Matrix<double,6,NV> F_t;
      
      // [ABA]
      typedef Eigen::Matrix<double,6,NV> U_t;
      typedef Eigen::Matrix<double,NV,NV> D_t;
      typedef Eigen::Matrix<double,6,NV> UD_t;

      typedef Eigen::Matrix<double,NQ,1> ConfigVector_t;
      typedef Eigen::Matrix<double,NV,1> TangentVector_t;
    };

  template<> struct traits<JointDataPrismaticUnaligned> { typedef JointPrismaticUnaligned JointDerived; };
  template<> struct traits<JointModelPrismaticUnaligned> { typedef JointPrismaticUnaligned JointDerived; };

  struct JointDataPrismaticUnaligned : public JointDataBase <JointDataPrismaticUnaligned>
  {
    typedef JointPrismaticUnaligned JointDerived;
    SE3_JOINT_TYPEDEF;

    Transformation_t M;
    Constraint_t S;
    Motion_t v;
    Bias_t c;

    F_t F;
    
    // [ABA] specific data
    U_t U;
    D_t Dinv;
    UD_t UDinv;

    JointDataPrismaticUnaligned() :
      M(1),
      S(Eigen::Vector3d::Constant(NAN)),
      v(Eigen::Vector3d::Constant(NAN),NAN),
      U(), Dinv(), UDinv()
    {}
    
    JointDataPrismaticUnaligned(const Motion_t::Vector3 & axis) :
      M(1),
      S(axis),
      v(axis,NAN),
      U(), Dinv(), UDinv()
    {}

  }; // struct JointDataPrismaticUnaligned

  struct JointModelPrismaticUnaligned : public JointModelBase <JointModelPrismaticUnaligned>
  {
    typedef JointPrismaticUnaligned JointDerived;
    SE3_JOINT_TYPEDEF;

    using JointModelBase<JointModelPrismaticUnaligned>::id;
    using JointModelBase<JointModelPrismaticUnaligned>::idx_q;
    using JointModelBase<JointModelPrismaticUnaligned>::idx_v;
    using JointModelBase<JointModelPrismaticUnaligned>::setIndexes;
    typedef Motion::Vector3 Vector3;
    
    JointModelPrismaticUnaligned() : axis(Vector3::Constant(NAN))   {}
    JointModelPrismaticUnaligned(Scalar x, Scalar y, Scalar z)
    {
      axis << x, y, z ;
      axis.normalize();
      assert(axis.isUnitary() && "Translation axis is not unitary");
    }
    
    JointModelPrismaticUnaligned(const Vector3 & axis) : axis(axis)
    {
      assert(axis.isUnitary() && "Translation axis is not unitary");
    }

    JointDataDerived createData() const { return JointDataDerived(axis); }
    
    void calc(JointDataDerived & data, const Eigen::VectorXd & qs) const
    {
      const double & q = qs[idx_q()];

      data.M.translation() = axis * q;
    }

    void calc(JointDataDerived & data,
              const Eigen::VectorXd & qs,
              const Eigen::VectorXd & vs) const
    {
      const Scalar & q = qs[idx_q()];
      const Scalar & v = vs[idx_v()];

      data.M.translation() = axis * q;
      data.v.rate = v;
    }
    
    void calc_aba(JointDataDerived & data, Inertia::Matrix6 & I, const bool update_I) const
    {
      data.U = I.block<6,3> (0,Inertia::LINEAR) * axis;
      data.Dinv[0] = 1./axis.dot(data.U.segment <3> (Inertia::LINEAR));
      data.UDinv.noalias() = data.U * data.Dinv;
      
      if (update_I)
        I -= data.UDinv * data.U.transpose();
    }
    
    Scalar finiteDifferenceIncrement() const
    {
      using std::sqrt;
      return sqrt(Eigen::NumTraits<Scalar>::epsilon());
    }

    ConfigVector_t neutralConfiguration_impl() const
    { 
      ConfigVector_t q;
      q << 0;
      return q;
    } 

    bool isSameConfiguration_impl(const Eigen::VectorXd& q1, const Eigen::VectorXd& q2, const Scalar & prec = Eigen::NumTraits<Scalar>::dummy_precision()) const
    {
      const Scalar & q_1 = q1[idx_q()];
      const Scalar & q_2 = q2[idx_q()];

      return (fabs(q_1 - q_2) < prec);
    }

    static std::string classname() { return std::string("JointModelPrismaticUnaligned"); }
    std::string shortname() const { return classname(); }

    Vector3 axis;
  }; // struct JointModelPrismaticUnaligned

} //namespace se3


#endif // ifndef __se3_joint_prismatic_unaligned_hpp__
