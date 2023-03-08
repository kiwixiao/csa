// STAR-CCM+ macro: ExportStarPlotData.java
// Written by STAR-CCM+ 14.06.012
package macro;

import java.util.*;

import star.common.*;
import star.base.neo.*;

public class ExportStarPlotData extends StarMacro {

private String subject = "Inspire04";
private String category = "CT003";

private String path = "/storage/Qiwei/cchmc_OSA/CCHMCProjects/PythonProjects/PlotsByPython";
private String outPath =path+"/files_"+subject+category+"/"+subject+category;


  public void execute() {
    execute0();
  }

  private void execute0() {

    Simulation simulation_0 = 
      getActiveSimulation();

    MonitorPlot monitorPlot_1 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("EnergyFlux_monitors_pts_Plot"));

    monitorPlot_1.export(resolvePath(outPath+"_EF.csv"), ",");

    MonitorPlot monitorPlot_2 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("EnergyFluxDiff_monitors_pts_Plot"));

    monitorPlot_2.export(resolvePath(outPath+"_EFDiff.csv"), ",");

    MonitorPlot monitorPlot_3 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("FdotV_Areamonitors_pts_Plot"));

    monitorPlot_3.export(resolvePath(outPath+"_AreaFdotV.csv"), ",");

    MonitorPlot monitorPlot_4 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("FdotV_Areamotionwise_monitors_pts_Plot"));

    monitorPlot_4.export(resolvePath(outPath+"_AreaFdotV_motionwise.csv"), ",");

    MonitorPlot monitorPlot_5 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("FdotV_AreamotionwiseAnterior_monitors_pts_Plot"));

    monitorPlot_5.export(resolvePath(outPath+"_AreaFdotV_motionwiseAnterior.csv"), ",");

    MonitorPlot monitorPlot_6 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("FdotV_AreamotionwiseLeft_monitors_pts_Plot"));

    monitorPlot_6.export(resolvePath(outPath+"_AreaFdotV_motionwiseLeft.csv"), ",");

    MonitorPlot monitorPlot_7 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("FdotV_AreamotionwisePosterior_monitors_pts_Plot"));

    monitorPlot_7.export(resolvePath(outPath+"_AreaFdotV_motionwisePosterior.csv"), ",");

    MonitorPlot monitorPlot_8 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("FdotV_AreamotionwiseRight_monitors_pts_Plot"));

    monitorPlot_8.export(resolvePath(outPath+"_AreaFdotV_motionwiseRight.csv"), ",");

    MonitorPlot monitorPlot_9 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("FdotV_monitors_pts_Plot"));

    monitorPlot_9.export(resolvePath(outPath+"_FdotV.csv"), ",");

    MonitorPlot monitorPlot_10 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("FdotV_motionwise_monitors_pts_Plot"));

    monitorPlot_10.export(resolvePath(outPath+"_FdotV_motionwise.csv"), ",");

    MonitorPlot monitorPlot_11 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("FdotV_motionwiseAnterior_monitors_pts_Plot"));

    monitorPlot_11.export(resolvePath(outPath+"_FdotV_motionwiseAnterior.csv"), ",");

    MonitorPlot monitorPlot_12 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("FdotV_motionwiseLeft_monitors_pts_Plot"));

    monitorPlot_12.export(resolvePath(outPath+"_FdotV_motionwiseLeft.csv"), ",");

    MonitorPlot monitorPlot_13 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("FdotV_motionwisePosterior_monitors_pts_Plot"));

    monitorPlot_13.export(resolvePath(outPath+"_FdotV_motionwisePosterior.csv"), ",");

    MonitorPlot monitorPlot_14 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("FdotV_motionwiseRight_monitors_pts_Plot"));

    monitorPlot_14.export(resolvePath(outPath+"_FdotV_motionwiseRight.csv"), ",");

    MonitorPlot monitorPlot_0 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("MassFlow_monitors_pts_Plot"));

    monitorPlot_0.export(resolvePath(outPath+"_MF.csv"), ",");

    MonitorPlot monitorPlot_15 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("Resistance_monitors_pts_Plot"));

    monitorPlot_15.export(resolvePath(outPath+"_Resis.csv"), ",");

    MonitorPlot monitorPlot_16 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("SurfaceArea_monitors_pts_Plot"));

    monitorPlot_16.export(resolvePath(outPath+"_SurfaceArea.csv"), ",");

    MonitorPlot monitorPlot_17 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("TotalPressure_monitors_pts_Plot"));

    monitorPlot_17.export(resolvePath(outPath+"_TP.csv"), ",");

    MonitorPlot monitorPlot_18 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("TotalPressureLoss_monitors_pts_Plot"));

    monitorPlot_18.export(resolvePath(outPath+"_TPL.csv"), ",");

    MonitorPlot monitorPlot_19 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("TotalVolume Monitor_pts Plot"));

    monitorPlot_19.export(resolvePath(outPath+"_TV.csv"), ",");

    MonitorPlot monitorPlot_20 = 
      ((MonitorPlot) simulation_0.getPlotManager().getPlot("VolumeSwept_monitors_pts_Plot"));

    monitorPlot_20.export(resolvePath(outPath+"_VS.csv"), ",");
  }
}
