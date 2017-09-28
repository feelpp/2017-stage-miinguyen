/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
*/

#include <feel/feel.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description bratuoptions( "Module Peltier Comsol" );
    bratuoptions.add_options()
    ( "k", Feel::po::value<double>()->default_value( 1.6 ), "_thermal_conductivity of SemiConductorP" )
    
    ( "alpha", Feel::po::value<double>()->default_value( 200e-6 ), "_seebeck_coeff of SemiConductorP" )
    
    ( "sigma", Feel::po::value<double>()->default_value( 1.1e5 ), "_electric_conductivity of SemiConductorP" )
    
    ( "V_N", Feel::po::value<double>()->default_value( -0.7e6 ), "Neumann condition of Intensity" ) //-0.7/(0.001*0.001)

    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
    
    ( "Peltier-Seebeck", Feel::po::value<bool>()->default_value( false ), "Ajoute effet Peltier - Seebeck   non=0, oui=1" )

    //auto matSemiConductorP = materialDesc(_seebeck_coeff=200e-6,_electric_conductivity=1.1e5,_thermal_conductivity=1.6,_mass_density=7740,_heat_capacity=154.4);
    //auto matElectrode = materialDesc(_seebeck_coeff=6.5e-6,_electric_conductivity=5.9e8,_thermal_conductivity=350,_mass_density=8920,_heat_capacity=385);
    
    ;
    return bratuoptions.add( Feel::feel_options() );
}

/**
 * Module Peltier Comsol
 *
 *
 */
int
main( int argc, char** argv )
{

    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="peltier_comsol",
                                  _author=" ",
                                  _email=" "));

    //typedef Mesh<Simplex<2>> mesh_type;
    typedef Mesh< Simplex< FEELPP_DIM,1 > > mesh_type;
    auto mesh = loadMesh(_mesh=new mesh_type);
    typedef FunctionSpace<mesh_type, bases<Lagrange<1,Scalar>,Lagrange<1,Scalar> > > space_type;
    auto Vh = space_type::New(_mesh=mesh);
    //auto Vh = THch< 1 >(mesh);
    auto TV = Vh->element();
    auto T = TV.element<0>();
    auto V = TV.element<1>();
    
    auto t = TV.element<0>();
    auto v = TV.element<1>();
    
    
    double k = doption(_name="k");
    
    double alpha = doption(_name="alpha");
    
    double sigma = doption(_name="sigma");
    
    double V_N = doption(_name="V_N");
    
    //auto V_N = expr( soption(_name="functions.V_N"), "V_N" ); //Neumann condition nonhomogene
    
    auto P = alpha*idv(T);
    
    auto d_P = alpha*idt(T);
    
    auto E = -gradv(V);
    auto d_E = -gradt(V);
    
    auto j1 = sigma*(E );
    auto d_j1 = sigma*(d_E );
    
    auto j2 = sigma*(E - alpha*gradv(T));
    auto d_j2 = sigma*(d_E - alpha*gradt(T));
    
    auto q1 = -k*gradv(T);
    auto q2 = -k*gradv(T) + P*j2;
    
    auto d_q1 = -k*gradt(T);
    auto d_q2 = -k*gradt(T)+(d_P*j2+P*d_j2);
    
    // Joule effect
    auto Q1 = sigma*inner( gradv(V) , gradv(V));
    auto Q2 = sigma*inner( gradv(V) + alpha*gradv(T) , gradv(V));
    
    auto d_Q1 = sigma*( inner( 2*gradt(V), gradv(V) ) );
    auto d_Q2 = sigma*( inner( 2*gradt(V) + alpha*gradt(T) , gradv(V) ) + alpha*inner( gradv(T) , gradt(V) ) );
    
    // electric charge equation
    auto d_fluxElectric1 = -sigma*(gradt(V) ); //d_j
    auto d_fluxElectric2 = -sigma*(gradt(V) + alpha*gradt(T));
    
    auto fluxElectric1 = -sigma*(gradv(V));
    auto fluxElectric2 = -sigma*(gradv(V)+alpha*gradv(T));
    

    auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
        {
            if (!J) J = backend()->newMatrix( Vh, Vh );
            auto l = form1(_test=Vh);
            auto a = form2( _test=Vh, _trial=Vh, _matrix=J );
            
            
            if ( boption("Peltier-Seebeck") )
            {
                // energy equation
                a = integrate( _range = elements(mesh), _expr= inner(d_q2,grad(t)) );
                
                a += integrate( _range = elements(mesh), _expr= -inner( d_Q2,id(t) ) );
                
                a +=on(_range=markedfaces(mesh,"Ground"),_rhs = l, _element=T, _expr = cst(0.) );
                
                // electric charge equation
                a += integrate( _range = elements(mesh),
                               _expr= inner(d_fluxElectric2 ,grad(v)) );
            
            }
            else
            {
                // energy equation
                a = integrate( _range = elements(mesh), _expr= inner(d_q1,grad(t)) );
                
                a += integrate( _range = elements(mesh), _expr= -inner( d_Q1,id(t) ) );
                
                a +=on(_range=markedfaces(mesh,"Ground"),_rhs = l, _element=T, _expr = cst(0.) );
                
                // electric charge equation
                a += integrate( _range = elements(mesh),
                               _expr= inner(d_fluxElectric1 ,grad(v)) );
            }
            
        };
    
    auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
        {
            auto TV = Vh->element();
            
            TV = *X;
            auto r = form1( _test=Vh, _vector=R );
            
            if ( boption("Peltier-Seebeck") )
            {
                // energy equation
                r = integrate( _range = elements(mesh), _expr= inner(q2,grad(t)) );
                r += integrate( _range = elements(mesh), _expr= -inner(Q2,id(T)) );
                // electric charge equation
                r += integrate(_range = elements(mesh), _expr= inner(fluxElectric2 ,grad(v)) );
            }
            else
            {
                // energy equation
                r = integrate( _range = elements(mesh), _expr= inner(q1,grad(t)) );
                r += integrate( _range = elements(mesh), _expr= -inner(Q1,id(t)) );
                // electric charge equation
                r += integrate( _range = elements(mesh), _expr= inner(fluxElectric1 ,grad(v)) );
            }
            
            r += integrate( _range=markedfaces(mesh,"Intensity"), _expr=  V_N * id(v)); //Neumann condition
            
            auto w = Vh->element();
            auto Tw = w.element<0>();
            auto Vw = w.element<1>();
            
            w=*R;
            
            Tw.on(_range=markedfaces(mesh,"Ground"),_expr = cst(0.) );
            Vw.on(_range=markedfaces(mesh,"Ground"),_expr = cst(0.) );

            *R=w;
            
        };
    T.on(_range = elements(mesh), _expr = cst(0.));
    V.on(_range = elements(mesh), _expr = cst(0.));
    //TV.zero();
    
    backend()->nlSolver()->residual = Residual;
    backend()->nlSolver()->jacobian = Jacobian;
    backend()->nlSolve( _solution= TV );

    auto e = exporter( _mesh=mesh );
    e->add( "T", T );
    e->add( "V", V );
    e->save();
}





