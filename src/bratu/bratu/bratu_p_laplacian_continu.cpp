/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */


#include <feel/feel.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description bratuoptions( "Bratu problem options" );
    bratuoptions.add_options()
    ( "p", Feel::po::value<double>()->default_value( 2 ), "p-type Laplacian's problem" )
    
    ( "step", Feel::po::value<double>()->default_value( 0.1 ), "step for increment reach to p" )
    
    ( "lambda", Feel::po::value<double>()->default_value( 1 ), "exp() coefficient value for the Bratu problem" )
    
    ( "epsilon", Feel::po::value<double>()->default_value( 0.00001 ), "coefficient of regulation" )
    
    ( "penalbc", Feel::po::value<double>()->default_value( 30 ), "penalisation parameter for the weak boundary conditions" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
    
    ( "export-matlab", "export matrix and vectors in matlab" )
    
    ( "weak-method", Feel::po::value<bool>()->default_value( false ), "weak/fort method for Dirichlet's boudary condition   fort=0, weak=1" )
    ;
    return bratuoptions.add( Feel::feel_options() );
}

/**
 * Bratu p-laplacian Problem
 *
 * solve \f$ −∇⋅(η∇u)−λexp(u) u = 0, \quad u_\Gamma = 0\f$ on \f$\Omega\f$
 avec $$\eta = (\epsilon^2+\gamma)^{\frac{p-2}{2}},\quad \gamma=\frac{1}{2}|\nabla u|^2$$
 */

int
main( int argc, char** argv )
{
    
    using namespace Feel;
    using Feel::cout;
    
    Environment env( _argc=argc, _argv=argv,
                    _desc=makeOptions(),
                    _about=about(_name="bratu_p-laplacian_3D",
                                 _author=" ",
                                 _email=" "));
    //auto mesh = unitSquare();
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    
    double p = doption(_name="p");
    double step = doption(_name="step");
    double penalbc = doption(_name="penalbc");
    double lambda = doption(_name="lambda");
    double epsilon = doption(_name="epsilon");
    
    u.zero();
    
    for (double i = 2; i <= (p+step); i+=step)
    {
        if (i != p)
        {
            auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
            {
                if (!J) J = backend()->newMatrix( Vh, Vh );
                auto a = form2( _test=Vh, _trial=Vh, _matrix=J );
                auto l = form1(_test=Vh);
                auto gamma = 0.5* inner(gradv(u), gradv(u));
                auto eta = pow( epsilon*epsilon + gamma, 0.5*(i - 2) );
                a = integrate( elements( mesh ), eta * (gradt( u )*trans( grad( v ) ) ));
                a += integrate( elements( mesh ), 0.5*(i - 2) * pow( epsilon*epsilon  + gamma, 0.5*(i - 4) ) *
                           ( inner(gradt(u),gradv(u)) ) *
                           ( gradv(u) * trans(grad(v)) )  );
                a += integrate( elements( mesh ), lambda*( exp( idv( u ) ) )*idt( u )*id( v ) );
                
                a +=on(_range=boundaryfaces( mesh ),_rhs = l, _element=u, _expr=cst(0.) );
            
            };
        
            auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
            {
                auto u = Vh->element();
                u = *X;
                auto r = form1( _test=Vh, _vector=R );
                auto gamma = 0.5* inner(gradv(u), gradv(u));
                auto eta = pow( epsilon*epsilon  + gamma, 0.5*(i - 2) );
                r = integrate( elements( mesh ), eta * ( gradv( u )*trans( grad( v ) ) ) );
                r += integrate( elements( mesh ),  lambda*exp( idv( u ) )*id( v ) );
                
                auto v = Vh->element();
                v=*R; // copy residual in v
                // set the unknowns on the boundary to 0
                v.on(_range=boundaryfaces(mesh),_expr=cst(0.));
                // copy back to R
                *R=v;
            
            };
        
            cout << "----------------------------------------------" <<std::endl;
            cout << " p = " << i <<std::endl;
            backend()->nlSolver()->residual = Residual;
            backend()->nlSolver()->jacobian = Jacobian;
            backend()->nlSolve( _solution=u );
            cout << "----------------------------------------------" <<std::endl;
        }
    }
    
    
    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();
}





