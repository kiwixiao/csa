// STAR-CCM+ macro: exportBoundaryForCSA.java
// Written by STAR-CCM+ 14.06.012
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;
import star.meshing.*;

public class exportBoundaryForCSA extends StarMacro {

  private String sub = "Inspire04";
  private String outPath="/storage/Qiwei/cchmc_OSA/CCHMCProjects/PythonProjects/CSAmeasurements/";
  private String category = "CT003";
  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();

    SurfaceRep surfaceRep_0 = 
      ((SurfaceRep) simulation_0.getRepresentationManager().getObject("Extracted Surface"));

    SurfaceRepRegion surfaceRepRegion_0 = 
      ((SurfaceRepRegion) surfaceRep_0.getSurfaceRepRegionManager().getObject("Region"));

    SurfaceRepBoundary Nasopharynx = 
      ((SurfaceRepBoundary) surfaceRepRegion_0.getSurfaceBoundaries().getObject(sub+category+".Nasopharynx"));
    SurfaceRepBoundary Oropharynx = 
      ((SurfaceRepBoundary) surfaceRepRegion_0.getSurfaceBoundaries().getObject(sub+category+".Oropharynx"));
    SurfaceRepBoundary Larynx = 
      ((SurfaceRepBoundary) surfaceRepRegion_0.getSurfaceBoundaries().getObject(sub+category+".Larynx"));
    SurfaceRepBoundary Trachea = 
      ((SurfaceRepBoundary) surfaceRepRegion_0.getSurfaceBoundaries().getObject(sub+category+".Trachea"));
    SurfaceRepBoundary LeftNose = 
      ((SurfaceRepBoundary) surfaceRepRegion_0.getSurfaceBoundaries().getObject(sub+category+".LeftNose"));
    SurfaceRepBoundary RightNose = 
      ((SurfaceRepBoundary) surfaceRepRegion_0.getSurfaceBoundaries().getObject(sub+category+".RightNose"));
    SurfaceRepBoundary Mask = 
      ((SurfaceRepBoundary) surfaceRepRegion_0.getSurfaceBoundaries().getObject(sub+category+".Mask"));

    Nasopharynx.exportStlBoundary(new NeoObjectVector(new Object[] {Nasopharynx}), resolvePath(outPath+sub+category+"/"+sub+category+"_Nasopharynx.stl"), true);

    Oropharynx.exportStlBoundary(new NeoObjectVector(new Object[] {Oropharynx}), resolvePath(outPath+sub+category+"/"+sub+category+"_Oropharynx.stl"), true);

    Larynx.exportStlBoundary(new NeoObjectVector(new Object[] {Larynx}), resolvePath(outPath+sub+category+"/"+sub+category+"_Larynx.stl"), true);

    Trachea.exportStlBoundary(new NeoObjectVector(new Object[] {Trachea}), resolvePath(outPath+sub+category+"/"+sub+category+"_Trachea.stl"), true);

    Mask.exportStlBoundary(new NeoObjectVector(new Object[] {Mask}), resolvePath(outPath+sub+category+"/"+sub+category+"_Mask.stl"), true);

    RightNose.exportStlBoundary(new NeoObjectVector(new Object[] {RightNose}), resolvePath(outPath+sub+category+"/"+sub+category+"_RightNose.stl"), true);

    LeftNose.exportStlBoundary(new NeoObjectVector(new Object[] {LeftNose,Nasopharynx,Oropharynx,Larynx,Trachea}), resolvePath(outPath+sub+category+"/"+sub+category+"_LeftNoseDecending.stl"), true);
    
  }
}
