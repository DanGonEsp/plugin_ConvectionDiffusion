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

template <int dim> struct face_type_traits
{
    typedef void face_type0;
	typedef void face_type1;
};

template <> struct face_type_traits<1>
{
    typedef ReferenceVertex face_type0;
	typedef ReferenceVertex face_type1;
};

template <> struct face_type_traits<2>
{
    typedef ReferenceEdge face_type0;
	typedef ReferenceEdge face_type1;
};

template <> struct face_type_traits<3>
{
    typedef ReferenceTriangle face_type0;
	typedef ReferenceQuadrilateral face_type1;
};

// parameters:
// secure_container& sides : as given by associated_elements collection computed beforehand
// const TGridFunction& u : grid function
// std::vector<MultiIndex<2> > multInd : constrained dof_indices computed beforehand
// size_t fct, specifies funtion used in multi index computation, default 0
template<typename side_type,typename secure_container,typename TGridFunction>
	void get_constrained_sides_cr(secure_container& sides,const TGridFunction& u,std::vector<MultiIndex<2> > multInd,size_t fct = 0){
	size_t nOfSides = sides.size();
	ConstrainingEdge* cEdge=NULL;
	ConstrainingFace* cFace=NULL;
	Edge* edge;
	Face* face;
	size_t nc;
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;
	///	world dimension
	static const int dim = domain_type::dim;
	///	create Multiindex
	std::vector<MultiIndex<2> > seMultInd;
	for (size_t i=0;i<nOfSides;i++){
		if (dim==2){
			cEdge = dynamic_cast<ConstrainingEdge*>(sides[i]);
			if (cEdge==NULL) continue;
			nc = cEdge->num_constrained_edges();
			cEdge->constrained_edge(0);
		}
		else{
			cFace = dynamic_cast<ConstrainingFace*>(sides[i]);
			if (cFace==NULL) continue;
			nc = cFace->num_constrained_faces();
		}
		for (size_t k=0;k<nc;k++){
			if (dim==2){
				edge = dynamic_cast<Edge*>(cEdge->constrained_edge(k));
				u.inner_dof_indices(edge,fct,seMultInd);
			} else {
				face = 	dynamic_cast<Face*>(cFace->constrained_face(k));
				u.inner_dof_indices(face,fct,seMultInd);
			}
			for(size_t j=nOfSides;j<multInd.size();j++){
				if (multInd[j][0]==seMultInd[0][0]){
					if (dim==2) sides.push_back(dynamic_cast<side_type*>(edge));
					else sides.push_back(dynamic_cast<side_type*>(face));
					break;
				}
			}
		}
	}
}

template <typename TGridFunction,typename side_type,typename constraining_side_type,typename VType>
void constrainingSideAveraging(PeriodicAttachmentAccessor<side_type,Attachment<VType> >& aaData,SmartPtr<TGridFunction> m_uInfo){
	//	domain type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::grid_type grid_type;
	static const int dim = domain_type::dim;
	/// side iterator
	typedef typename TGridFunction::template traits<constraining_side_type>::const_iterator cSideIterator;
	typedef typename domain_type::position_accessor_type position_accessor_type;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;
	domain_type& domain = *m_uInfo->domain().get();
	DimCRFVGeometry<dim> geo;
	position_accessor_type posAcc = m_uInfo->domain()->position_accessor();
	cSideIterator cSideIter = m_uInfo->template begin<constraining_side_type>(SurfaceView::SHADOW_RIM);
	cSideIterator cSideIterEnd = m_uInfo->template end<constraining_side_type>(SurfaceView::SHADOW_RIM);
	for(  ;cSideIter !=cSideIterEnd; ++cSideIter){
		constraining_side_type* cSide = *cSideIter;
		typename grid_type::template traits<elem_type>::secure_container assoElements;
		typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
		side_secure_container sides;
		// get associated element
		domain.grid()->associated_elements(assoElements, cSide);
		elem_type* elem = assoElements[0];
		domain.grid()->associated_elements_sorted(sides, elem);
		std::vector<DoFIndex> ind;
		m_uInfo->dof_indices(elem,0,ind,true,true);
		get_constrained_sides_cr<side_type,side_secure_container,TGridFunction>(sides,*m_uInfo,ind);
		elem = assoElements[0];
		std::vector<MathVector<dim> > vCorner;
		CollectCornerCoordinates(vCorner, *elem, posAcc);
		geo.update_hanging(elem, &(vCorner[0]), domain.subset_handler().get());
		for (size_t i=0;i<geo.num_constrained_dofs();i++){
			const typename DimCRFVGeometry<dim>::CONSTRAINED_DOF& cd = geo.constrained_dof(i);
			const size_t index = cd.index();
			if (dynamic_cast<side_type*>(sides[index])!=dynamic_cast<side_type*>(*cSideIter)) continue;
			aaData[sides[index]]*=0;
			for (size_t j=0;j<cd.num_constraining_dofs();j++){
				VType localValue;
				size_t cdIndex = cd.constraining_dofs_index(j);
				localValue = aaData[sides[cdIndex]];
				localValue *= cd.constraining_dofs_weight(j);
				aaData[sides[index]] += localValue;
			}
		}
	}
}


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

	/// blockSize of used algebra
		static const int blockSize = algebra_type::blockSize;

	///	grid type
		typedef typename domain_type::grid_type grid_type;

	/// element type
		typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// side type
		typedef typename elem_type::side side_type;

	/// element iterator
		typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// side iterator
		typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;

		static const size_t _P_ = dim;
        static const size_t _C_ = dim+1;

	///	Type of geometric base object
		typedef typename domain_traits<TDomain::dim>::grid_base_object grid_base_object;

	/// position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
	
		typedef std::vector<std::pair<DoFIndex, MathVector<dim> > > vIndexPosPair;
		typedef std::vector<std::vector<std::pair<DoFIndex, MathVector<dim> > > > vvIndexPosPair;
		typedef std::pair<MathVector<dim>, MathVector<dim> > MathVector_Pair;
		
		typedef MathMatrix<dim,dim> dimMat;
		typedef Attachment<dimMat> AMathDimMat;
		typedef PeriodicAttachmentAccessor<side_type,AMathDimMat > aSideDimMat;
		typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
	
		typedef Attachment<std::vector< MathVector<dim> > > ANumberArray;
		typedef Attachment<std::vector< DoFIndex > > ASizetArray;
		typedef PeriodicAttachmentAccessor<side_type,ANumberArray> aSideNumberArray;
		typedef PeriodicAttachmentAccessor<side_type,ASizetArray> aSideSizetArray;
	
		aSideDimMat acGrad;
		aSideNumber acVol;
		aSideNumberArray acGradSh;
		aSideSizetArray acGradShInd;
		
		AMathDimMat aGrad;
		ANumber aVol;
		ANumberArray aGradSh;
		ASizetArray aGradShInd;
		
		typedef typename face_type_traits<dim>::face_type0 face_type0;
		typedef typename face_type_traits<dim>::face_type1 face_type1;

	private:
		SmartPtr<TGridFunction> m_u;
		grid_type* m_grid;
		bool m_bAdaptive;
		bool m_bLinPressureDefect;
		bool m_bLinPressureJacobian;
		bool m_bLinUpConvDefect;
		bool m_bLinUpConvJacobian;
		bool m_limiter;
        number m_BackFLowValue;
		ISubsetHandler* m_ish;
		// zero gradient subset group
		SubsetGroup m_zeroGradSg;
		
	public:
		void init(SmartPtr<TGridFunction> u, number vBackFlowValue){
			m_u = u;
            m_BackFLowValue=vBackFlowValue;
			domain_type& domain = *m_u->domain().get();
			grid_type& grid = *domain.grid();
			m_grid = &grid;
			m_ish = m_u->domain()->subset_handler().get();

		}
		


	/// constructor
        TrasportConstraintFVC(SmartPtr<TGridFunction> u, number vBackFlowValue){
			init(u,vBackFlowValue);
		};


	///	destructor
		~TrasportConstraintFVC() {};



	
		// add linear pressure part and linear velocity upwind part to jacobian
		virtual void adjust_jacobian(matrix_type& J, const vector_type& u,
				                             ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
				                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = NULL,const number s_a0 = 1.0){
            /*
			domain_type& domain = *m_u->domain().get();

			//	create Multiindex
			std::vector<DoFIndex> multInd;
			
			position_accessor_type aaPos = m_u->domain()->position_accessor();
			
			typename grid_type::template traits<side_type>::secure_container sides;
			std::vector<MathVector<dim> > vCorner;
			std::vector<DoFIndex> ind;
			
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
					
					/// handle convection
					if (true){
					
						MathVector<dim> vVel;

                        for(size_t ip = 0; ip < geo.num_scv();; ++ip)
						{
							const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(ip);
                            size_t s = scv.node_id();
                            
							VecSet(vVel, 0.0);
                            for (int d=0;d<dim;d++){
                                dd->inner_dof_indices(sides[s], d, ind);//Possible error
                                vVel[d] = DoFRef(u,ind[0]);
                            }
                            number prod = s_a0*VecDot(vVel,scv.normal());
                            
                            
                            m_grid->associated_elements(assoElements,sides[s]);
                            size_t numOfAsso = assoElements.size();
                            
                            if(numOfAsso==1)
                            {
                                if(prod>0){
                                    
                                    dd->dof_indices(elem, _C_ , multInd);
                                    DoFRef(d,multInd[0]) += flux;
                                }
                            }
                            else//numOfAsso=2
                            {
                                if(prod>0)
                                {
                                    if (assoElements[0]==elem)
                                    {
                                        dd->dof_indices(assoElements[0], _C_ , multInd);
                                        DoFRef(d,multInd[0]) += flux;
                                        DoFRef(J,localInd,shapeInd)+= flux;
                                        
                                        dd->dof_indices(assoElements[1], _C_ , multInd);
                                        DoFRef(d,multInd[0]) -= flux;
                                    }
                                    else
                                    {
                                        dd->dof_indices(assoElements[0], _C_ , multInd);
                                        DoFRef(d,multInd[0]) -= flux;
                                        
                                        dd->dof_indices(assoElements[1], _C_ , multInd);
                                        DoFRef(d,multInd[0]) += flux;
                                        
                                        
                                        DoFRef(J,localInd,shapeInd)+= flux;
                                        localInd[0]=ind[scvf.to()][0]+d0;
                                        DoFRef(J,localInd,shapeInd)-= flux;
                                    }
                                }
                            }
                            
                            
                            
                            

						}

						for(size_t ip = 0; ip < geo.num_scv(); ++ip){
	

							side_type* baseSide = sides[base];
							MathVector<dim> distVec;
							VecSubtract(distVec, scvf.global_ip(),geo.scv(base).global_ip());
							for	(int d0=0;d0<dim;d0++){
								size_t shapeSize = acGradShInd[baseSide].size();
								for (size_t sh=0;sh<shapeSize;sh++){
									for (int d1=0;d1<dim;d1++){
										number flux = prod * distVec[d1] * acGradSh[baseSide][sh][d1];
										localInd[0]=ind[scvf.from()][0]+d0;
										shapeInd[0]=acGradShInd[baseSide][sh][0]+d0;
										DoFRef(J,localInd,shapeInd)+= flux;
										localInd[0]=ind[scvf.to()][0]+d0;
										DoFRef(J,localInd,shapeInd)-= flux;
									}
								}	
							}
						} 
					}// if m_bLinUpConvJacobian
					
					/// handle pressure
					if (m_bLinPressureJacobian==true){
					
						//  get sides of element
						m_grid->associated_elements_sorted(sides, elem );
					
						//	reference object type
						ReferenceObjectID roid = elem->reference_object_id();
					
						//	compute size (volume) of element
						const number elemSize = ElementSize<dim>(roid, &vCorner[0]);
					
						typename grid_type::template traits<elem_type>::secure_container assoElements;
					
						std::vector<DoFIndex> elemInd(sides.size()+1);
					
						dd->inner_dof_indices(elem,_P_,ind);
						elemInd[0] = ind[0];
					
						//UG_LOG("0 " << elemInd[0] << "\n");

						size_t gradShapesSize = 1;
					
						//UG_LOG(elem << "\n");

						MathMatrix<2*dim+1,dim> gradShapes;
						for (int sh=0;sh<2*dim+1;sh++) for (int d0=0;d0<dim;d0++) gradShapes[sh][d0]=0;

						for (size_t s=0;s<sides.size();s++){
							m_grid->associated_elements(assoElements,sides[s]);
							// face value is average of associated elements
							size_t numOfAsso = assoElements.size();
							const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(s);
							if (numOfAsso==1){
								for (int d=0;d<dim;d++) gradShapes[0][d]+=scv.normal()[d];
								continue;
							}
							for (size_t i=0;i<numOfAsso;i++){
								dd->inner_dof_indices(assoElements[i],_P_,ind);
								//UG_LOG(assoElements[i] << "\n");
								if (assoElements[i]!=elem) break;
							}
							elemInd[gradShapesSize] = ind[0];
							//UG_LOG(gradShapesSize << " " << elemInd[gradShapesSize] << "\n");
							for (int d=0;d<dim;d++){
								gradShapes[0][d]+=0.5*scv.normal()[d];
								gradShapes[gradShapesSize][d]+=0.5*scv.normal()[d];
							}
							gradShapesSize++;
						}
						gradShapes/=(number)elemSize;
						//UG_LOG("gradShapesSize=" << gradShapesSize << "\n");
						//UG_LOG(gradShapes << "\n");
						// for debug set grad shapes to trivial
						//for (int d2=0;d2<dim;d2++) gradShapes[0][d2]=1;
						//for (size_t sh=0;sh<gradShapesSize;sh++)for (int d2=0;d2<dim;d2++) gradShapes[sh][d2]=0;
					

						for (int d1=0;d1<dim;d1++){
							dd->dof_indices(elem, d1 , multInd);
							for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
							{
								// 	get current SCVF
								const typename DimCRFVGeometry<dim>::SCVF& scvf = geo.scvf(ip);
								MathVector<dim> distVec;
								VecSubtract(distVec,scvf.global_ip(),geo.global_bary());

								for (size_t sh=0;sh<gradShapesSize;sh++){
									for (int d2=0;d2<dim;d2++){
										number flux = s_a0 * distVec[d2]*gradShapes[sh][d2];
										//UG_LOG("from = " << multInd[scvf.from()] << " p = " << elemInd[sh] << "\n");
										DoFRef(J,multInd[scvf.from()],elemInd[sh])+= flux *  scvf.normal()[d1];
										//UG_LOG("to = " << multInd[scvf.to()] << " p = " << elemInd[sh] << "\n");
										DoFRef(J,multInd[scvf.to()],elemInd[sh])-= flux *  scvf.normal()[d1];
									}
								}
							}
						}
					}//if m_bLinPressureJacobian
				}//For Iter
			}*/
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
                add_convection_defect(d,u,dd);
            }
			else {
				//	loop all time points and assemble them
				for(size_t t = 0; t < vScaleStiff->size(); ++t){
					if ((*vScaleStiff)[t]==0) continue;
                    else
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
