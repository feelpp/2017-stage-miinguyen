/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
*/

#include <feel/feel.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description bratuoptions( "Module Peltier Comsol" );
    bratuoptions.add_options()
    ( "k1", Feel::po::value<double>()->default_value( 1.6 ), "_thermal_conductivity of SemiConductorP" )
    ( "k2", Feel::po::value<double>()->default_value( 350 ), "_thermal_conductivity of Electrode" )
    
    ( "alpha1", Feel::po::value<double>()->default_value( 200e-6 ), "_seebeck_coeff of SemiConductorP" )
    ( "alpha2", Feel::po::value<double>()->default_value( 6.5e-6 ), "_seebeck_coeff of Electrode" )
    
    ( "sigma1", Feel::po::value<double>()->default_value( 1.1e5 ), "_electric_conductivity of SemiConductorP" )
    ( "sigma2", Feel::po::value<double>()->default_value( 5.9e8 ), "_electric_conductivity of Electrode" )
    
    ( "V_N", Feel::po::value<double>()->default_value( -0.7e6 ), "Neumann condition of Intensity" ) //-0.7/(0.001*0.001)

    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )

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
    using Feel::cout;
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
    
    
    double k1 = doption(_name="k1");
    double k2 = doption(_name="k2");
    
    double alpha1 = doption(_name="alpha1");
    double alpha2 = doption(_name="alpha2");
    
    double sigma1 = doption(_name="sigma1");
    double sigma2 = doption(_name="sigma2");
    
    double V_N = doption(_name="V_N");
    
    //auto V_N = expr( soption(_name="functions.V_N"), "V_N" ); //Neumann condition nonhomogene
    
    auto P1 = alpha1*idv(T);
    auto P2 = alpha2*idv(T);
    
    auto d_P1 = alpha1*idt(T);
    auto d_P2 = alpha2*idt(T);
    
    auto E = -gradv(V);
    auto d_E = -gradt(V);
    
    auto j1 = sigma1*(E - alpha1*gradv(T));
    auto d_j1 = sigma1*(d_E - alpha1*gradt(T));
    
    auto j2 = sigma2*(E - alpha2*gradv(T));
    auto d_j2 = sigma2*(d_E - alpha2*gradt(T));
    
    auto q1 = -k1*gradv(T) + P1*j1;
    auto q2 = -k2*gradv(T) + P2*j2;
    
    auto d_q1 = -k1*gradt(T)+(d_P1*j1+P1*d_j1);
    auto d_q2 = -k2*gradt(T)+(d_P2*j2+P2*d_j2);
    
    // Joule effect
    auto Q1 = sigma1*inner( gradv(V) + alpha1*gradv(T) , gradv(V));
    auto Q2 = sigma2*inner( gradv(V) + alpha2*gradv(T) , gradv(V));
    
    auto d_Q1 = sigma1*( inner( 2*gradt(V) + alpha1*gradt(T) , gradv(V) ) + alpha1*inner( gradv(T) , gradt(V) ) );
    auto d_Q2 = sigma2*( inner( 2*gradt(V) + alpha2*gradt(T) , gradv(V) ) + alpha2*inner( gradv(T) , gradt(V) ) );
    
    // electric charge equation
    auto d_fluxElectric1 = -sigma1*(gradt(V) + alpha1*gradt(T)); //d_j
    auto d_fluxElectric2 = -sigma2*(gradt(V) + alpha2*gradt(T));
    
    auto fluxElectric1 = -sigma1*(gradv(V)+alpha1*gradv(T));
    auto fluxElectric2 = -sigma2*(gradv(V)+alpha2*gradv(T));
    

    auto Jacobian = [&](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
        {
            if (!J) J = backend()->newMatrix( Vh, Vh );
            auto l = form1(_test=Vh);
            auto a = form2( _test=Vh, _trial=Vh, _matrix=J );
            
            TV = *X;
            
            // energy equation
            a = integrate( _range=markedelements(mesh, "Material0"), _expr= -inner(d_q1,grad(t)) );
            a += integrate( _range=markedelements(mesh, "Material1"), _expr= -inner(d_q2,grad(t)) );
            
            a += integrate( _range=markedelements(mesh,"Material0"), _expr= -inner( d_Q1,id(t) ) );
            a += integrate( _range=markedelements(mesh,"Material1"), _expr= -inner( d_Q2,id(t) ) );
            
            
            // electric charge equation
            a += integrate( _range=markedelements(mesh,"Material0"),
                           _expr= -inner(d_fluxElectric1 ,grad(v)) ); //
            a += integrate( _range=markedelements(mesh,"Material1"),
                           _expr= -inner(d_fluxElectric2 ,grad(v)) ); //
            
            a +=on(_range=markedfaces(mesh,"Ground"),_rhs = l, _element=T, _expr = cst(0.) );
            a += on(_range=markedfaces(mesh,"Ground"),_rhs = l, _element=V, _expr = cst(0.) );
            
        };
    
    auto Residual = [&](const vector_ptrtype& X, vector_ptrtype& R)
        {
            //auto TV = Vh->element();
            
            TV = *X;
            auto r = form1( _test=Vh, _vector=R );
            
            // energy equation
            r = integrate( _range=markedelements(mesh,"Material0"), _expr= -inner(q1,grad(t)) );
            r += integrate( _range=markedelements(mesh,"Material1"), _expr= -inner(q2,grad(t)) );
            
            r += integrate( _range=markedelements(mesh,"Material0"), _expr= -inner(Q1,id(t)) );
            r += integrate( _range=markedelements(mesh,"Material1"), _expr= -inner(Q2,id(t)) );
            
            // electric charge equation
            r += integrate( _range=markedelements(mesh,"Material0"), _expr= -inner(fluxElectric1 ,grad(v)) );
            r += integrate( _range=markedelements(mesh,"Material1"), _expr= -inner(fluxElectric2 ,grad(v)) );
            
            r += integrate( _range=markedfaces(mesh,"Intensity"), _expr=  V_N * id(v)); //Neumann condition
            
            auto w = Vh->element();
            auto Tw = w.element<0>();
            auto Vw = w.element<1>();
            
            w=*R;
            
            Tw.on(_range=markedfaces(mesh,"Ground"),_expr = cst(0.) );
            Vw.on(_range=markedfaces(mesh,"Ground"),_expr = cst(0.) );

            *R=w;
            
        };
    T.on(_range = elements(mesh), _expr = cst(273.15));
    V.on(_range = elements(mesh), _expr = cst(0.));
    
    backend()->nlSolver()->residual = Residual;
    backend()->nlSolver()->jacobian = Jacobian;
    backend()->nlSolve( _solution= TV );
    
    auto T1 = mean(_range = markedelements(mesh,"Electrode1"), _expr  = idv(T))(0,0);
    
    cout << "T_electrode1 = " << T1<< std::endl;
    
    auto T2 = mean(_range = markedelements(mesh,"Electrode2"), _expr  = idv(T))(0,0);
    
    cout << "T_electrode2 = " << T2<< std::endl;
    
    auto T3 = mean(_range = markedelements(mesh,"Material1"), _expr  = idv(T))(0,0);
    
    cout << "T_Adiabatic = " << T3<< std::endl;
    
    auto T4 = mean(_range = markedelements(mesh), _expr  = idv(T))(0,0) ;
    
    cout << "T_mean = " << T2<< std::endl;
    
    
    // Exporter le gradient de V
    auto XhVec = Pdhv<1>(mesh);
    auto electricField = XhVec->element();
    electricField.on(_range=elements(mesh),_expr=trans(gradv(V)));
    
    auto e = exporter( _mesh=mesh );
    e->add( "T", T );
    e->add( "V", V );
    e->add( "electric-field", electricField );
    e->save();
}





