/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__PLUGINS__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FVC__NEIGHBOR__CONCENTRATION__
#define __H__UG__PLUGINS__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FVC__NEIGHBOR__CONCENTRATION__

#include "common/common.h"

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton_update_interface.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_grid/tools/subset_group.h"
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "lib_grid/algorithms/attachment_util.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace ConvectionDiffusionPlugin{




/**
concept derived from grid_function_user_data.h
 */
template <typename TGridFunction>
class UpwindValueSource
:     public StdUserData<UpwindValueSource<TGridFunction>,
      number, TGridFunction::dim>,
      virtual public INewtonUpdate
      {
    ///    domain type
    typedef typename TGridFunction::domain_type domain_type;

    ///    algebra type
    typedef typename TGridFunction::algebra_type algebra_type;

    /// position accessor type
    typedef typename domain_type::position_accessor_type position_accessor_type;

    ///    world dimension
    static const int dim = domain_type::dim;

    ///    grid type
    typedef typename domain_type::grid_type grid_type;

    /// element type
    typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

      /// side type
      typedef typename elem_type::side side_type;
    /// MathVector<dim> attachment
    //        typedef MathVector<dim> vecDim;
    //        typedef Attachment<vecDim> AMathVectorDim;

    /// attachment accessor
    typedef PeriodicAttachmentAccessor<side_type,ANumber > aSideNumber;
    typedef Grid::AttachmentAccessor<elem_type,ANumber > aElementNumber;

    /// element iterator
    typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

    /// Side iterator
    typedef typename TGridFunction::template traits<side_type>::const_iterator SideIterator;
    
    typedef  DimCRFVGeometry<dim> TFVGeo;
          
    typedef typename grid_type::template traits<side_type>::secure_container side_secure_container;
    typedef typename grid_type::template traits<elem_type>::secure_container elem_secure_container;

          private:
    // old pressure attachment accessor
    ANumber m_aPOld;
    aElementNumber m_pOld;

    //    pressure attachment accessor (interpolated pressure in Sides)
    ANumber m_aP;
    aSideNumber m_p;

    // level set grid function
    SmartPtr<TGridFunction> m_u;

    //    approximation space for level and surface grid
    SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

    //  grid
    grid_type* m_grid;

    //  pressure index
    size_t _C_;

          private:

    ///    Data import for source
    SmartPtr<CplUserData<MathVector<dim>,dim> > m_imSource;

          public:
    /////////// Source

    void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
    {
        m_imSource = data;
    }

    void set_source(number f_x)
    {
        SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
        for (int i=0;i<dim;i++){
            f->set_entry(i, f_x);
        }
        set_source(f);
    }

    void set_source(number f_x, number f_y)
    {
        if (dim!=2){
            UG_THROW("NavierStokes: Setting source vector of dimension 2"
                    " to a Discretization for world dim " << dim);
        } else {
            SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
            f->set_entry(0, f_x);
            f->set_entry(1, f_y);
            set_source(f);
        }
    }

    void set_source(number f_x, number f_y, number f_z)
    {
        if (dim<3){
            UG_THROW("NavierStokes: Setting source vector of dimension 3"
                    " to a Discretization for world dim " << dim);
        }
        else
        {
            SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
            f->set_entry(0, f_x);
            f->set_entry(1, f_y);
            f->set_entry(2, f_z);
            set_source(f);
        }
    }

#ifdef UG_FOR_LUA
    void set_source(const char* fctName)
    {
        set_source(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
    }
#endif

          public:
    /// constructor
        UpwindValueSource(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct, size_t DOF){
        
        _C_=DOF;
        m_u = spGridFct;
        domain_type& domain = *m_u->domain().get();
        grid_type& grid = *domain.grid();
        m_grid = &grid;
        m_spApproxSpace = approxSpace;
        set_source(0.0);
        grid.template attach_to<elem_type>(m_aPOld);
        grid.template attach_to<side_type>(m_aP);
        m_pOld.access(grid,m_aPOld);
        m_p.access(grid,m_aP);
      
        // set all values to zero
        SetAttachmentValues(m_p, m_u->template begin<side_type>(), m_u->template end<side_type>(), 0);
        SetAttachmentValues(m_pOld, m_u->template begin<elem_type>(), m_u->template end<elem_type>(), 0);
        this->update();
    }

    virtual ~UpwindValueSource(){};

    template <int refDim>
    inline void evaluate(number vValue[],
                         const MathVector<dim> vGlobIP[],
                         number time, int si,
                         GridObject* elem,
                         const MathVector<dim> vCornerCoords[],
                         const MathVector<refDim> vLocIP[],
                         const size_t nip,
                         LocalVector* u,
                         const MathMatrix<refDim, dim>* vJT = NULL) const
    {

        std::vector<DoFIndex> multInd;
        
        side_secure_container sides;
        
        UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
        
        elem_type* element = static_cast<elem_type*>(elem);
        
        m_grid->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );
        


        const size_t numVertices = element->num_vertices();
        //    get domain of grid function
        const domain_type& domain = *m_u->domain().get();

        //    get position accessor

        const position_accessor_type& posAcc = domain.position_accessor();

//        position_accessor_type aaPos = m_u->domain()->position_accessor();

        // coord and Side pos array
        MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
        Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
        DimCRFVGeometry<dim> crfvgeo;

        for(size_t i = 0; i < numVertices; ++i){
            vVrt[i] = element->vertex(i);
            coCoord[i] = posAcc[vVrt[i]];
        };
        

        // evaluate finite volume geometry
        crfvgeo.update(elem, &(coCoord[0]), domain.subset_handler().get());
        
        size_t nSides = crfvgeo.num_scv();

        std::vector<MathVector<dim>> vVel(nSides);
        

        
        for (size_t i=0;i < nSides;i++)
        {
            const typename DimCRFVGeometry<dim>::SCV& scv = crfvgeo.scv(i);
            size_t s = scv.node_id();

            for (int d=0;d<dim;d++){
                m_u->dof_indices(sides[i], d, multInd);
                vVel[i][d] = DoFRef(*m_u,multInd[0]);
                
            }
            if(VecDot(vVel[i],scv.normal())>0)
            {
                m_u->inner_dof_indices(element, _C_, multInd);
                vValue[s]=DoFRef(*m_u,multInd[0]);
            }
            else
            {
                elem_secure_container elems;
                m_grid->associated_elements(elems, static_cast<side_type*>(sides[i]));
                if(elems.size()==1)
                {
                    vValue[s]=0;
                }
                else
                {
                    if (elems[0]==element)
                    {
                        m_u->inner_dof_indices(elems[1], _C_, multInd);
                        vValue[s]=DoFRef(*m_u,multInd[0]);
                    }
                    else
                    {
                        m_u->inner_dof_indices(elems[0], _C_, multInd);
                        vValue[s]=DoFRef(*m_u,multInd[0]);
                    }
                }
            }
        }

    }; // evaluate

    void update(){
        //    get domain
        domain_type& domain = *m_u->domain().get();

        //    create Multiindex
        std::vector<DoFIndex> multInd;
        
        DimCRFVGeometry<dim> geo;


        SetAttachmentValues(m_p, m_u->template begin<side_type>(), m_u->template end<side_type>(), 0);

        
        
        //    coord and vertex array
        MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
        Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

        //    get position accessor
        const position_accessor_type& posAcc = domain.position_accessor();
        
        for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
            //    get iterators
            ElemIterator iter = m_u->template begin<elem_type>(si);
            ElemIterator iterEnd = m_u->template end<elem_type>(si);

            //    loop elements of dimension
            for(  ;iter !=iterEnd; ++iter)
            {
                //    get Elem
                elem_type* elem = *iter;
                
                //    get vertices and extract corner coordinates
                const size_t numVertices = elem->num_vertices();
                for(size_t i = 0; i < numVertices; ++i){
                    vVrt[i] = elem->vertex(i);
                    coCoord[i] = posAcc[vVrt[i]];
                };

                static const size_t MaxNumSidesOfElem = 10;

                typedef MathVector<dim> MVD;
                std::vector<MVD> uValue(MaxNumSidesOfElem);
 
                side_secure_container sides;

                UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Only elements of type elem_type are currently supported");

                domain.grid()->associated_elements_sorted(sides, static_cast<elem_type*>(elem) );

                //    evaluate finite volume geometry

                geo.update(elem, &(coCoord[0]), domain.subset_handler().get());


                size_t nNaturalSides = geo.num_sh();
                size_t nSides = geo.num_scv();

                for (size_t s=0;s < nNaturalSides;s++){
                    for (int d=0;d<dim;d++){
                        m_u->dof_indices(sides[s], d, multInd);
                        uValue[s][d]=DoFRef(*m_u,multInd[0]);
                    }
                }
                m_u->inner_dof_indices(elem, _C_, multInd);
                m_pOld[elem]=DoFRef(*m_u,multInd[0]);
                
                
    
                for (size_t i=0;i < nSides;i++)
                {
                    const typename DimCRFVGeometry<dim>::SCV& scv = geo.scv(i);
                    size_t s = scv.node_id();
                    if(VecDot(uValue[s],scv.normal())>0)
                    {
                        m_p[sides[s]] += m_pOld[elem];
                    }
                
                }

            }
        }
    }

          private:
    static const size_t max_number_of_ips = 20;

          public:
    virtual void operator() (number& value,
                             const MathVector<dim>& globIP,
                             number time, int si) const
    {
        UG_THROW("UpwindValue: Need element.");
    }

    virtual void operator() (number vValue[],
                             const MathVector<dim> vGlobIP[],
                             number time, int si, const size_t nip) const
    {
        UG_THROW("LevelSetUserData: Need element.");
    }

          
          
    virtual void compute(LocalVector* u, GridObject* elem,
                         const MathVector<dim> vCornerCoords[], bool bDeriv = false)
    {
        const int si = this->subset();
        for(size_t s = 0; s < this->num_series(); ++s)
            evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
                          elem, NULL, this->template local_ips<dim>(s),
                          this->num_ip(s), u);
    }

    virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
                         const MathVector<dim> vCornerCoords[], bool bDeriv = false)
    {
        const int si = this->subset();
        for(size_t s = 0; s < this->num_series(); ++s)
            evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
                          elem, NULL, this->template local_ips<dim>(s),
                          this->num_ip(s), &(u->solution(this->time_point(s))));
    }

    ///    returns if provided data is continuous over geometric object boundaries
    virtual bool continuous() const {return false;}

    ///    returns if grid function is needed for evaluation
    virtual bool requires_grid_fct() const {return true;}
};




/**
concept derived from grid_function_user_data.h
 */
template <typename TGridFunction>
class GradientSource
: 	public StdUserData<GradientSource<TGridFunction>, MathVector<TGridFunction::dim>, TGridFunction::dim>,
  	virtual public INewtonUpdate
  	{
	///	domain type
	typedef typename TGridFunction::domain_type domain_type;

	///	algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	///	world dimension
	static const int dim = domain_type::dim;

	///	grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// MathVector<dim> attachment
	//		typedef MathVector<dim> vecDim;
	//		typedef Attachment<vecDim> AMathVectorDim;

	/// attachment accessor
	typedef PeriodicAttachmentAccessor<Vertex,ANumber > aVertexNumber;
	typedef Grid::AttachmentAccessor<elem_type,ANumber > aElementNumber;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

  		private:
	// old pressure attachment accessor
	ANumber m_aPOld;
	aElementNumber m_pOld;

	//	pressure attachment accessor (interpolated pressure in vertices)
	ANumber m_aP;
	aVertexNumber m_p;

	//  volume attachment accessor
	ANumber m_aVol;
	aVertexNumber m_vol;

	// level set grid function
	SmartPtr<TGridFunction> m_u;

	//	approximation space for level and surface grid
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

	//  grid
	grid_type* m_grid;

	//  pressure index
	size_t _C_;

  		private:

	///	Data import for source
	SmartPtr<CplUserData<MathVector<dim>,dim> > m_imSource;

  		public:
	/////////// Source

	void set_source(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
	{
		m_imSource = data;
	}

	void set_source(number f_x)
	{
		SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
		for (int i=0;i<dim;i++){
			f->set_entry(i, f_x);
		}
		set_source(f);
	}

	void set_source(number f_x, number f_y)
	{
		if (dim!=2){
			UG_THROW("NavierStokes: Setting source vector of dimension 2"
					" to a Discretization for world dim " << dim);
		} else {
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			set_source(f);
		}
	}

	void set_source(number f_x, number f_y, number f_z)
	{
		if (dim<3){
			UG_THROW("NavierStokes: Setting source vector of dimension 3"
					" to a Discretization for world dim " << dim);
		}
		else
		{
			SmartPtr<ConstUserVector<dim> > f(new ConstUserVector<dim>());
			f->set_entry(0, f_x);
			f->set_entry(1, f_y);
			f->set_entry(2, f_z);
			set_source(f);
		}
	}

#ifdef UG_FOR_LUA
	void set_source(const char* fctName)
	{
		set_source(LuaUserDataFactory<MathVector<dim>, dim>::create(fctName));
	}
#endif

  		public:
	/// constructor
        GradientSource(SmartPtr<ApproximationSpace<domain_type> > approxSpace,SmartPtr<TGridFunction> spGridFct, size_t DOF){
        
        _C_=DOF;
		m_u = spGridFct;
		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		m_spApproxSpace = approxSpace;
		set_source(0.0);
		grid.template attach_to<elem_type>(m_aPOld);
		grid.template attach_to<Vertex>(m_aP);
		grid.template attach_to<Vertex>(m_aVol);
		m_pOld.access(grid,m_aPOld);
		m_p.access(grid,m_aP);
		m_vol.access(grid,m_aVol);
		// set all values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_p, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_pOld, m_u->template begin<elem_type>(), m_u->template end<elem_type>(), 0);
		this->update();
	}

	virtual ~GradientSource(){};

	template <int refDim>
	inline void evaluate(MathVector<dim> vValue[],
	                     const MathVector<dim> vGlobIP[],
	                     number time, int si,
	                     GridObject* elem,
	                     const MathVector<dim> vCornerCoords[],
	                     const MathVector<refDim> vLocIP[],
	                     const size_t nip,
	                     LocalVector* u,
	                     const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);

		//	reference object id
		ReferenceObjectID roid = elem->reference_object_id();

		const size_t numVertices = element->num_vertices();
		//    get domain of grid function
		const domain_type& domain = *m_u->domain().get();

		//    get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

//		position_accessor_type aaPos = m_u->domain()->position_accessor();

		// coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];
		DimCRFVGeometry<dim> crfvgeo;

		MathVector<dim> grad[domain_traits<dim>::MaxNumVerticesOfElem];

		for(size_t i = 0; i < numVertices; ++i){
			vVrt[i] = element->vertex(i);
			coCoord[i] = posAcc[vVrt[i]];
		};

		for (size_t i=0;i<domain_traits<dim>::MaxNumVerticesOfElem;i++)
			grad[i]*=0.0;

		// evaluate finite volume geometry
		crfvgeo.update(elem, &(coCoord[0]), domain.subset_handler().get());

		// Lagrange 1 trial space
		const LocalShapeFunctionSet<dim>& lagrange1 =
				LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

		std::vector<number> shapes;

		for (size_t ip=0;ip<crfvgeo.num_scvf();ip++){
			const typename DimCRFVGeometry<dim>::SCVF& scvf = crfvgeo.scvf(ip);
			number pinter=0;
			lagrange1.shapes(shapes,scvf.local_ip());
			for (size_t sh=0;sh<numVertices;sh++)
				pinter += m_p[vVrt[sh]]*shapes[sh];
			for (int d=0;d<dim;d++){
				number flux = pinter*scvf.normal()[d];
				grad[scvf.from()][d]-=flux;
				grad[scvf.to()][d]+=flux;
			}
		}

		for (size_t i=0;i<crfvgeo.num_scv();i++){
			grad[i]/=crfvgeo.scv(i).volume();
		}

		// evaluate source data
		(*m_imSource)(vValue,
				vGlobIP,
				time, si,
				elem,
				vCornerCoords,
				vLocIP,
				nip,
				u,
				vJT);
		for (size_t i=0;i<nip;i++)
			vValue[i]+=grad[i];
	}; // evaluate

	void update(){
		//	get domain
		domain_type& domain = *m_u->domain().get();
		//	create Multiindex
		std::vector<DoFIndex> multInd;
		DimFV1Geometry<dim> geo;

		//	coord and vertex array
		MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];
		Vertex* vVrt[domain_traits<dim>::MaxNumVerticesOfElem];

		//	get position accessor
		typedef typename domain_type::position_accessor_type position_accessor_type;
		const position_accessor_type& posAcc = domain.position_accessor();

		// set volume and p values to zero
		SetAttachmentValues(m_vol, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		SetAttachmentValues(m_p, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0);
		// set p^{old} = p^{old} + p^{sep}
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			ElemIterator iter = m_u->template begin<elem_type>(si);
			ElemIterator iterEnd = m_u->template end<elem_type>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				elem_type* elem = *iter;
				m_u->inner_dof_indices(elem, _C_, multInd);
				m_pOld[elem]=DoFRef(*m_u,multInd[0]);
			}
		}
		// compute pressure in vertices by averaging
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			ElemIterator iter = m_u->template begin<elem_type>(si);
			ElemIterator iterEnd = m_u->template end<elem_type>(si);
			for(  ;iter !=iterEnd; ++iter)//Loop for estimate pressure at all vertices
            {
				elem_type* elem = *iter;
				number pValue = m_pOld[elem];
				const size_t numVertices = elem->num_vertices();
				for(size_t i = 0; i < numVertices; ++i){
					vVrt[i] = elem->vertex(i);
					coCoord[i] = posAcc[vVrt[i]];
				};
				geo.update(elem, &(coCoord[0]), domain.subset_handler().get());
				for(size_t i = 0; i < numVertices; ++i){
					number scvVol = geo.scv(i).volume();
					m_vol[vVrt[i]]+=scvVol;
					m_p[vVrt[i]]+=scvVol*pValue;
				}
                
			}
		}
		PeriodicBoundaryManager* pbm = (domain.grid())->periodic_boundary_manager();
		// go over all vertices and average
		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si){
			VertexIterator iter = m_u->template begin<Vertex>(si);
			VertexIterator iterEnd = m_u->template end<Vertex>(si);
			for(  ;iter !=iterEnd; ++iter)
			{
				Vertex* vrt = *iter;
				if (pbm && pbm->is_slave(vrt)) continue;
				m_p[vrt]/=m_vol[vrt];
			}
		}
	}

  		private:
	static const size_t max_number_of_ips = 20;

  		public:
	virtual void operator() (MathVector<dim>& value,
	                         const MathVector<dim>& globIP,
	                         number time, int si) const
	{
		UG_THROW("LevelSetUserData: Need element.");
	}

	virtual void operator() (MathVector<dim> vValue[],
	                         const MathVector<dim> vGlobIP[],
	                         number time, int si, const size_t nip) const
	{
		UG_THROW("LevelSetUserData: Need element.");
	}

	virtual void compute(LocalVector* u, GridObject* elem,
	                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
			              elem, NULL, this->template local_ips<dim>(s),
			              this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
	                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si,
			              elem, NULL, this->template local_ips<dim>(s),
			              this->num_ip(s), &(u->solution(this->time_point(s))));
	}

	///	returns if provided data is continuous over geometric object boundaries
	virtual bool continuous() const {return false;}

	///	returns if grid function is needed for evaluation
	virtual bool requires_grid_fct() const {return true;}
};


} // namespace CONVENCTIONDIFFUSION
} // end namespace ug


#endif /*__H__UG__PLUGINS__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FVC__NEIGHBOR__CONCENTRATION__*/
