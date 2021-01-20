/*
 * RandomWall.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: fionnoh <fionnoh@gmail.com>
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

#ifndef Hadrons_MSource_RandomWall_hpp_
#define Hadrons_MSource_RandomWall_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Random Wall source
 -----------------------------
 
 * options:
 - tW: source timeslice (integer)
 
 */

/******************************************************************************
 *                         Random Wall                                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class RandomWallPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RandomWallPar,
                                    unsigned int, tW);
};

template <typename FImpl>
class TRandomWall: public Module<RandomWallPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TRandomWall(const std::string name);
    // destructor
    virtual ~TRandomWall(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool        hasT_{false};
    std::string tName_;
};

// MODULE_REGISTER_TMP(RandomWall, TRandomWall<FIMPL>, MSource);
MODULE_REGISTER_TMP(StagRandomWall, TRandomWall<STAGIMPL>, MSource);

/******************************************************************************
 *                 TRandomWall implementation                                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TRandomWall<FImpl>::TRandomWall(const std::string name)
: Module<RandomWallPar>(name)
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TRandomWall<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TRandomWall<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandomWall<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envTmpLat(LatticeComplex, "eta1");
    envTmpLat(LatticeComplex, "eta2");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandomWall<FImpl>::execute(void)
{    
    LOG(Message) << "Generating random wall source at t = " << par().tW 
                 << std::endl;
    
    auto  &src = envGet(PropagatorField, getName());
    auto  &t   = envGet(Lattice<iScalar<vInteger>>, tName_);
    auto  nc   = FImpl::Dimension;
    
    if (!hasT_)
    {
        LatticeCoordinate(t, Tp);
        hasT_ = true;
    }

    envGetTmp(LatticeComplex, eta1);
    envGetTmp(LatticeComplex, eta2);
    src = Zero();
    for (unsigned int c = 0; c < nc; ++c)
    {
        gaussian(rng4d(), eta1);
        gaussian(rng4d(), eta2);
        eta1 += timesI(eta2);
        eta1 = where((t == par().tW), eta1, 0.*eta1);
        pokeColour(src, eta1, c,c);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_RandomWall_hpp_
