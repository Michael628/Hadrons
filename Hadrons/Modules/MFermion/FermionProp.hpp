/*
 * FermionProp.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Guido Cossu <guido.cossu@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Nils Asmussen <n.asmussen@soton.ac.uk>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: pretidav <david.preti@csic.es>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */

#ifndef Hadrons_MFermion_FermionProp_hpp_
#define Hadrons_MFermion_FermionProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                FermionProp                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class FermionPropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FermionPropPar,
                                    std::string, source,
                                    std::string, solver);
};

template <typename FImpl>
class TStagFermionProp: public Module<FermionPropPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TStagFermionProp(const std::string name);
    // destructor
    virtual ~TStagFermionProp(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
    Solver       *solver_{nullptr};
};

MODULE_REGISTER_TMP(StagFermionProp, TStagFermionProp<STAGIMPL>, MFermion);

/******************************************************************************
 *                      TStagFermionProp implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStagFermionProp<FImpl>::TStagFermionProp(const std::string name)
: Module<FermionPropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStagFermionProp<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().solver};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TStagFermionProp<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_5d"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagFermionProp<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().solver);
    envCreateLat(FermionField, getName());
    envTmpLat(FermionField, "tmp");
    if (Ls_ > 1)
    {
        envTmpLat(FermionField, "src", Ls_);
        envTmpLat(FermionField, "sol", Ls_);
        envCreateLat(FermionField, getName() + "_5d", Ls_);
    }
    else
    {
        envTmpLat(FermionField, "src");
        envTmpLat(FermionField, "sol");
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagFermionProp<FImpl>::execute(void)
{
    LOG(Message) << "Computing quark propagator '" << getName() << "'"
    << std::endl;
    
    std::string propName = (Ls_ == 1) ? getName() : (getName() + "_5d");
    auto        &prop    = envGet(FermionField, propName);
    auto        &source = envGet(FermionField, par().source);
    auto        &solver  = envGet(Solver, par().solver);
    auto        &mat     = solver.getFMat();
    
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, src);
    envGetTmp(FermionField, tmp);
    LOG(Message) << "Inverting using solver '" << par().solver
    << "' on source '" << par().source << "'" << std::endl;

    LOG(Message) << "Import source" << std::endl;
    if (!env().isObject5d(par().source))
    {
        if (Ls_ == 1) {
            src = source;
        } else {
            tmp = source;
            mat.ImportPhysicalFermionSource(tmp,src);
        }
    }
    // source conversion for 5D sources
    else
    {
        if (Ls_ != env().getObjectLs(par().source))
        {
            HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
        } else {
            src = source;
        }
    }
    LOG(Message) << "Solve" << std::endl;
    //sol = zero;
    sol = Zero();
    solver(prop, src);
    LOG(Message) << "Export solution" << std::endl;
    //std::cout<< "color " << c << " sol= " << sol << std::endl;
    // create 4D propagators from 5D one if necessary
    if (Ls_ > 1)
    {
        FermionField &p4d = envGet(FermionField, getName());
        mat.ExportPhysicalFermionSolution(prop, p4d);
    }
}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_FermionProp_hpp_
