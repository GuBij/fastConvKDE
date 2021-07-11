//#include "Uforest.h"
#include "UopenField.h"
#include "KTaylorOF.h"
#include "UhomTurb.h"
#include "KTaylorHT.h"
#include "UCouette.h"
#include "KTaylorCF.h"
//#include "Kforest.h"
//#include "KopenField.h"
//#include "KTaylorForest.h"

 meteorology::meteorology ()
 {}

 meteorology* meteorology::New(bool LangevinModel, bool master_proc )
 {
   dictionary dict({"meteofile","terrainModel","stackHeight","diameter","flowRate","gasTemp","zMax","measureFreq","humidityType","latitude","windVarType","dtByTauL","minTimeStep"});

   if ( dict.get(1) == "openField" )
	return new openField(dict, LangevinModel, master_proc);
   else if ( dict.get(1) == "CouetteFlow" )
	return new CouetteFlow(dict, LangevinModel, master_proc);
//   else if ( dict.get(1) == "forest" )
//        return new Uforest(dict, LangevinModel);
   else
	return nullptr;
 }

 meteorology* meteorology::New(string type, bool LangevinModel, bool master_proc )
 {
   dictionary dict({"meteofile","terrainModel","stackHeight","diameter","flowRate","gasTemp","zMax","measureFreq","humidityType","latitude","windVarType","dtByTauL","minTimeStep"});

   string modelName = dict.get(1);
   if (type == "windModel")
   {
	if (modelName == "openField")
	   return new UopenField(dict, master_proc);
//	else if (modelName == "forest")
//	   return new Uforest(dict);
	else if (modelName == "homTurb")
	   return new UhomTurb(dict, master_proc);
	else if (modelName == "CouetteFlow")
	   return new UCouette(dict, master_proc);
	else
	   return nullptr;
   }
   else if (type == "eddyDiffModel")
   {
        if (modelName == "openField" && LangevinModel)
	   return new KTaylorOF(dict, master_proc);
	else if (modelName == "homTurb" && LangevinModel)
	   return new KTaylorHT(dict, master_proc);
	else if (modelName == "CouetteFlow" && LangevinModel)
	   return new KTaylorCF(dict, master_proc);
/*
	else if (modelName == "openField")
           return new KopenField(dict);
	else if (modelName == "forest" && LangevinModel)
	   return new KTaylorForest(dict);
        else if (modelName == "forest")
           return new Kforest(dict);
*/
	else
	   return nullptr;
   }
   else
	error msg("inputDict",type);
 }

