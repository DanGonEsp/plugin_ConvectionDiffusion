/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
 * Author: Christian Wehner
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__PLUGINS__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FVC__TRANSPORT_CONSTRAINT_H_
#define __H__UG__PLUGINS__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FVC__TRANSPORT_CONSTRAINT_H_

#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"

namespace ug{


template <typename TGridFunction>
class TrasportConstraintFVC: public IDomainConstraint<typename TGridFunction::domain_type, typename TGridFunction::algebra_type>
{
	public:
		typedef typename TGridFunction::domain_type TDomain;
		typedef typename TGridFunction::algebra_type TAlgebra;

	///	world Dimension
		static const int dim = TDomain::dim;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	///	Type of Domain
		typedef TDomain domain_type;

	///	grid type
		typedef typename domain_type::grid_type grid_type;

	/// element type
		typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// side type
		typedef typename elem_type::side side_type;

	/// element iterator
		typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

        static const size_t _C_ = dim+1;

	/// position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;

	private:
		SmartPtr<TGridFunction> m_u;
		grid_type* m_grid;
		bool m_bFVC_ConvectionDefect;
		bool m_bFVC_Convection_lin_Defect;
    
        number m_BackFLowValue;

		
	public:
		void init(SmartPtr<TGridFunction> u, bool bConvectionDefect, bool bConvection_lin_Defect, number vBackFlowValue){
			m_u = u;
            m_bFVC_ConvectionDefect=bConvectionDefect;
            m_bFVC_Convection_lin_Defect=bConvection_lin_Defect;
            m_BackFLowValue=vBackFlowValue;
			domain_type& domain = *m_u->domain().get();
			grid_type& grid = *domain.grid();
			m_grid = &grid;
		}
		

	/// constructor
        TrasportConstraintFVC(SmartPtr<TGridFunction> u, bool bConvectionDefect, bool bConvection_lin_Defect, number vBackFlowValue){
			init(u, bConvectionDefect, bConvection_lin_Defect, vBackFlowValue);
		};


	///	destructor
		~TrasportConstraintFVC() {};



	
		// add linear pressure part and linear velocity upwind part to jacobian
		virtual void adjust_jacobian(matrix_type& J, const vector_type& u,
				                             ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
				                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = NULL,const number s_a0 = 1.0){
            
            if(!m_bFVC_ConvectionDefect) return;
            
			domain_type& domain = *m_u->domain().get();
			position_accessor_type aaPos = m_u->domain()->position_accessor();
			
			typename grid_type::template traits<side_type>::secure_container sides;
            
			std::vector<MathVector<dim> > vCorner;
            
            //    create Multiindex
            std::vector<DoFIndex> multInd;
			std::vector<DoFIndex> ind;
            std::vector<DoFIndex> ind2;
			
            DoFIndex localInd;
            localInd[1]=0;
            DoFIndex shapeInd;
            shapeInd[1]=0;
			//	create a FV Geometry for the dimension
			DimCRFVGeometry<dim> geo;

			for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
			{
			//	get iterators
				ElemIterator iter = dd->template begin<elem_type>(si);
				ElemIterator iterEnd = dd->template end<elem_type>(si);

			//	loop elements of dimension
				for(  ;iter !=iterEnd; ++iter)
				{
					//	get Elem
					elem_type* elem = *iter;
					
					//  get sides of element
					m_grid->associated_elements_sorted(sides, elem );
					
					//	get corners of element
					CollectCornerCoordinates(vCorner, *elem, aaPos);
					
					//	evaluate finite volume geometry
					geo.update(elem, &(vCorner[0]), domain.subset_handler().get());
                    
                    typename grid_type::template traits<elem_type>::secure_container assoElements;
					
					/// handle convection
					if (m_bFVC_ConvectionDefect){
					
						MathVector<dim> vVel;
                        MathVector<dim> vLin_Def;
                        
                        dd->inner_dof_indices(elem,_C_,ind);
                        number elemValue = DoFRef(u,ind[0]);

                        for(size_t ip = 0; ip < geo.num_scv(); ++ip)
						{
							const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(ip);
                            size_t s = scv.node_id();
                            
							VecSet(vVel, 0.0);
                            VecSet(vLin_Def, 0.0);
                            for (int d=0;d<dim;d++){
                                dd->dof_indices(elem,0,ind2);
                                localInd[0]=ind2[s][0]+d;
                                vVel[d] = DoFRef(u,localInd);
                            }
                            number flux = s_a0*VecDot(vVel,scv.normal());
                            VecScale(vLin_Def,scv.normal(),elemValue*s_a0);
                            
                            m_grid->associated_elements(assoElements,sides[s]);
                            size_t numOfAsso = assoElements.size();
                            
                            if(numOfAsso==1)
                            {
                                if(flux>0){
                                    
                                    dd->dof_indices(elem, _C_ , multInd);
                                    DoFRef(J,multInd[0],multInd[0]) += flux;
                                    
                                    if (m_bFVC_Convection_lin_Defect)
                                        for (int d=0;d<dim;d++){
                                            localInd[0]=ind2[s][0]+d;
                                            DoFRef(J,multInd[0],localInd) += vLin_Def[d];
                                        }
                                    
                                }
                                else
                                {
                                    if (m_bFVC_Convection_lin_Defect)
                                    {
                                        dd->dof_indices(elem, _C_ , multInd);
                                        for (int d=0;d<dim;d++){
                                            localInd[0]=ind2[s][0]+d;
                                            DoFRef(J,multInd[0],localInd) -= s_a0*m_BackFLowValue*scv.normal()[d];
                                        }
                                    }
                                }
                            }
                            else//numOfAsso=2
                            {
                                if(flux>0)
                                {
                                    if (assoElements[0]==elem)
                                    {
                                        dd->dof_indices(assoElements[0], _C_ , multInd);
                                        DoFRef(J,multInd[0],multInd[0]) += flux;
                                    
                                        dd->dof_indices(assoElements[1], _C_ , ind);
                                        DoFRef(J,ind[0],multInd[0]) -= flux;
                                        
                                        if (m_bFVC_Convection_lin_Defect)
                                        {
                                            
                                            for (int d=0;d<dim;d++){
                                                localInd[0]=ind2[s][0]+d;
                                                
                                                dd->dof_indices(assoElements[0], _C_ , multInd);
                                                DoFRef(J,multInd[0],localInd) += vLin_Def[d];
                                                
                                                dd->dof_indices(assoElements[1], _C_ , multInd);
                                                DoFRef(J,multInd[0],localInd) -= vLin_Def[d];
                                            }
                                        }
                                    }
                                    else
                                    {
                                        dd->dof_indices(assoElements[1], _C_ , multInd);
                                        DoFRef(J,multInd[0],multInd[0]) += flux;
                                        
                                        dd->dof_indices(assoElements[0], _C_ , ind);
                                        DoFRef(J,ind[0],multInd[0]) -= flux;
                                        
                                        if (m_bFVC_Convection_lin_Defect)
                                        {
                                            for (int d=0;d<dim;d++){
                                                localInd[0]=ind2[s][0]+d;
                                                
                                                dd->dof_indices(assoElements[1], _C_ , multInd);
                                                DoFRef(J,multInd[0],localInd) += vLin_Def[d];
                                                
                                                dd->dof_indices(assoElements[0], _C_ , multInd);
                                                DoFRef(J,multInd[0],localInd) -= vLin_Def[d];
                                                
                                            }
                                        }
                                    }
                                }
                            }
						}

					}// if m_bFVC_ConvectionDefect
					
				}//For Iter
			}
		};
		
        // add linear pressure part to defect
        virtual void add_convection_defect(vector_type& d, const vector_type& u, ConstSmartPtr<DoFDistribution>dd, const number time = 0.0,const number s_a = 1.0){
            
            domain_type& domain = *m_u->domain().get();

            //    create Multiindex
            std::vector<DoFIndex> multInd;
            
            position_accessor_type aaPos =  m_u->domain()->position_accessor();
            
            typedef typename grid_type::template traits<side_type>::secure_container secure_container;
            
            secure_container sides;
            std::vector<MathVector<dim> > vCorner;
            std::vector<DoFIndex> ind;

            //    create a FV Geometry for the dimension
            DimCRFVGeometry<dim> geo;

            for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
            {
                ElemIterator iter = dd->template begin<elem_type>(si);
                ElemIterator iterEnd = dd->template end<elem_type>(si);
                
            //    loop elements of dimension
                for(  ;iter !=iterEnd; ++iter)
                {
                    //    get Elem
                    elem_type* elem = *iter;

                    //  get sides of element
                    m_grid->associated_elements_sorted(sides, elem );
                    
                    //    get corners of element
                    CollectCornerCoordinates(vCorner, *elem, aaPos);
                    
                    //    evaluate finite volume geometry
                    geo.update(elem, &(vCorner[0]), domain.subset_handler().get());
                    
                    
                    typename grid_type::template traits<elem_type>::secure_container assoElements;
                    
                    dd->inner_dof_indices(elem,_C_,ind);
                    number elemValue = DoFRef(u,ind[0]);
                    
                    MathVector<dim> vVel;
                    
                    for (size_t i=0;i<geo.num_scv();i++)
                    {
                        
                        VecSet(vVel, 0);
                        const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(i);
                        size_t s = scv.node_id();
                        
                        for (int d=0;d<dim;d++){
                            dd->inner_dof_indices(sides[s], d, ind);//Possible error
                            vVel[d] = DoFRef(u,ind[0]);
                        }
                        number prod = s_a*VecDot(vVel,scv.normal());
                        number flux=elemValue*prod;
                        
                        m_grid->associated_elements(assoElements,sides[s]);
                        size_t numOfAsso = assoElements.size();
                        
                        if(numOfAsso==1)
                        {
                            if(prod<0)
                                flux=m_BackFLowValue*prod;
                            
                            dd->dof_indices(elem, _C_ , multInd);
                            DoFRef(d,multInd[0]) += flux;
                        }
                        else//numOfAsso=2
                        {
                            if(prod>0)
                            {
                                if (assoElements[0]==elem)
                                {
                                    dd->dof_indices(assoElements[0], _C_ , multInd);
                                    DoFRef(d,multInd[0]) += flux;
                                    
                                    dd->dof_indices(assoElements[1], _C_ , multInd);
                                    DoFRef(d,multInd[0]) -= flux;
                                }
                                else
                                {
                                    dd->dof_indices(assoElements[0], _C_ , multInd);
                                    DoFRef(d,multInd[0]) -= flux;
                                    
                                    dd->dof_indices(assoElements[1], _C_ , multInd);
                                    DoFRef(d,multInd[0]) += flux;
                                }
                            }
                        }
                    }//ForSCV
                }
            }
        };

			///	adapts defect to enforce constraints
			/// \{
		virtual void adjust_defect(vector_type& d, const vector_type& u,
				                           ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
				                           ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
										   const std::vector<number>* vScaleMass = NULL,
										   const std::vector<number>* vScaleStiff = NULL)
		{
            if (vSol == SPNULL){
                if(m_bFVC_ConvectionDefect)
                    add_convection_defect(d,u,dd);
            }
			else {
				//	loop all time points and assemble them
				for(size_t t = 0; t < vScaleStiff->size(); ++t){
					if ((*vScaleStiff)[t]==0) continue;
                    else
                        if(m_bFVC_ConvectionDefect)
                            add_convection_defect(d,*(vSol->solution(t)),dd,time,(*vScaleStiff)[t]);
        
				}
			}
		};
			/// \}

			///	adapts matrix and rhs (linear case) to enforce constraints
			/// \{
		virtual void adjust_linear(matrix_type& mat, vector_type& rhs,
				                           ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0){};
			/// \}

			///	adapts a rhs to enforce constraints
			/// \{
		virtual void adjust_rhs(vector_type& rhs, const vector_type& u,
				                        ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0){};
			/// \}

			///	sets the constraints in a solution vector
			/// \{
		virtual void adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd,
				                             int type, number time = 0.0){};

	///	returns the type of the constraints 
		virtual int type() const {return CT_CONSTRAINTS;}
};

} // end namespace ug


#endif /* __H__UG__PLUGINS__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FVC__TRANSPORT_CONSTRAINT_H_ */
