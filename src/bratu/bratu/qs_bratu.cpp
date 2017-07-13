/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:set fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4
*/


/*inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description bratuoptions( "Bratu problem options" );
    bratuoptions.add_options()
    ( "lambda", Feel::po::value<double>()->default_value( 1 ),
     "exp() coefficient value for the Bratu problem" )
    ( "penalbc", Feel::po::value<double>()->default_value( 30 ),
     "penalisation parameter for the weak boundary conditions" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ),
     "first h value to start convergence" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return bratuoptions.add( Feel::feel_options() );
}*/


 /**
 * Bratu Problem
 *
 * solve \f$ -\Delta u + \lambda \exp(u) = 0, \quad u_\Gamma = 0\f$ on \f$\Omega\f$
 */
#include <feel/feel.hpp>

int main( int argc, char** argv )
{
    
    using namespace Feel;
    using Feel::cout;
    po::options_description bratuoptions( "Bratu problem options" );
    bratuoptions.add_options()
    ( "stabilisation", po::value<bool>()->default_value( false ), "enable/disable stabilisation disable=0, enable=1" )
    ( "no-solve", po::value<bool>()->default_value( false ), "No solve" )
    ( "lambda", Feel::po::value<double>()->default_value( 1 ),
     "exp() coefficient value for the Bratu problem" )
    ;
    
    Environment env( _argc=argc, _argv=argv,
                    _desc=bratuoptions,
                    _about=about(_name="bratu",
                                 _author="Christophe Prud'homme",
                                 _email="christophe.prudhomme@feelpp.org"));
    
    tic();
    //auto mesh = unitSquare();
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<FEELPP_DIM,1>>);
    toc("loadMesh");
    
    tic();
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    
    //double penalbc = option(_name="penalbc").as<double>();
    double lambda = doption(_name="lambda");
    toc("Vh");
    
    tic();
    auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
    {
        auto a = form2( _test=Vh, _trial=Vh, _matrix=J );
        a = integrate( elements( mesh ), gradt( u )*trans( grad( v ) ) );
        a += integrate( elements( mesh ), lambda*( exp( idv( u ) ) )*idt( u )*id( v ) );
        a += integrate(_range=markedfaces(mesh,"Dirichlet"), _expr=cst(0.) );
    };
    toc("Jacobien");
    
    tic();
    auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
    {
        auto u = Vh->element();
        u = *X;
        auto r = form1( _test=Vh, _vector=R );
        r = integrate( elements( mesh ), gradv( u )*trans( grad( v ) ) );
        r +=  integrate( elements( mesh ),  lambda*exp( idv( u ) )*id( v ) );
    };
    toc("Residual");
    
    tic();
    u.zero();
    backend()->nlSolver()->residual = Residual;
    backend()->nlSolver()->jacobian = Jacobian;
    backend()->nlSolve( _solution=u );
    toc("solve");
                       
    auto e = exporter( _mesh=mesh );
    e->add( "uh", u );
    e->save();
    return 0;
}
