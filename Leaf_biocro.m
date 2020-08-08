function LeafA=Leaf_biocro(BioCro_params,Einput,Eio)

WeatherRH=BioCro_params('rh');
WeatherTemperature=BioCro_params('temp');
Air_CO2=BioCro_params('Catm');
WeatherWind=BioCro_params('windspeed');
Radiation_PAR=BioCro_params('incident_par_micromol');
Radiation_NIR=BioCro_params('radiation_nir');
Radiation_LW=BioCro_params('radiation_longwave');
PhotosynthesisType=BioCro_params('photosynthesis_type');
Vcmax25=BioCro_params('Vmax');
Jmax25=BioCro_params('Jmax');
Pressure=BioCro_params('atmospheric_pressure');
PhotosynthesisTheta=BioCro_params('theta');
Rd25=BioCro_params('Rd');
BallBerryIntercept=BioCro_params('b0');
BallBerrySlope=BioCro_params('b1');
Air_O2=BioCro_params('O2');
WaterStressFactor=BioCro_params('StomataWS');

%%% from EnvInput.txt and/or LeafMetaRun.m %%%
global pcfactor;  
ProteinTotalRatio=0.973;
pcfactor=1/ProteinTotalRatio;

GRNC=1.0;
%%% end %%%

PhotosynthesQ10=0;
R=8.314472E-3;%Gas constant KJ mole^{-1} K^{-1}
% Convert=1E6/(2.35E5); %Convert W m^{-2} to u moles m^{-2} s^{-1}
Convert=1; %BioCro passes in umoles m^{-2} s^{-1}
Boltzman=5.6697E-8; % Stefan-Boltzmann constant W m^{-2} K^{-4}
LatentHeatVaporization=44000.0;%J mole^{-1}
% Pressure=101325.0; % Standard atmospheric pressure Pa
ConstantsCp=29.3;
% PhotosynthesisTheta=0.76;
% Rd25=0.6; % WY 201810
% BallBerryIntercept=0.008;
% BallBerrySlope=10.5;
%Air_CO2=400.0;
% Air_O2=210.0;
% % WaterStressFunction=3;
% WaterStressFactor=1.0;
    
MaxError = 0.25; MinError = 0.001;
ErrorCount = 0; MaxErrorCount = 1;
Previous2Gs = 0; %Stomatal conductance moles/m2 leaf area/s
Relax = 0.0; % Relaxation value for oscillation
%Initialize variables
PreviousLeafState_Ci = 0;
Previous2Gs = 0; %Oscillation
PreviousLeafState_Gs = 0;
PreviousLeafState_Temperature =0;
PreviousLeafMassFlux_GrossAssimilation = 0;
PreviousLeafMassFlux_NetAssimilation = 0;
ErrorLeafState_Ci=1;
ErrorLeafState_Gs=1;  
ErrorLeafState_Temperature=1;
LeafTemperature= WeatherTemperature;% Initial Leaf temperature C
Ci = 0.7 * Air_CO2;%Initial Ci u moles/mole
Gs = 0.01; % Initial stomatal conductance moles/m2 leaf area/s
Gb = 10.2; % Initial boundary layer conductance moles/m2 leaf area/s

  % Rearrange Einput to use the same order as sample
  inE = importdata('MeM_input5_0.txt');
  GeneOrder = inE.textdata(2:size(inE.textdata, 1), 1);
  GeneNo = size(GeneOrder, 1);
  GeneMap = containers.Map(Einput{:, 1}, Einput{:, 2});
  ExpValue = ones(GeneNo, 1);
  for i=1:GeneNo
    if isKey(GeneMap, GeneOrder{i, 1})
      ExpValue(i, 1) = GeneMap(GeneOrder{i, 1});
    else
      fprintf("Gene %s is missing.\n", GeneOrder{i, 1});
    end
  end

    %Convergence loop for leaf
	while ((abs(ErrorLeafState_Ci) >= MinError || abs(ErrorLeafState_Gs) >= MinError ||abs(ErrorLeafState_Temperature) >= MinError) && ErrorCount <= MaxErrorCount)

	  PhotosynthesisRate=ComputPhotosynthesisRate(PhotosynthesisType,PhotosynthesQ10,Vcmax25,Jmax25,Rd25,R,LeafTemperature,Convert, Radiation_PAR,PhotosynthesisTheta,Ci,Air_O2,GRNC,ExpValue,Eio{:,2});
            NetAssimilation=PhotosynthesisRate(1);
            GrossAssimilation=PhotosynthesisRate(2);
            Rd=PhotosynthesisRate(3);
            GammaStar=PhotosynthesisRate(4);
        BoundaryCound=ComputeBoundaryLayerConductance(WeatherTemperature,LeafTemperature,WeatherRH,WeatherWind,Pressure,Gs,NetAssimilation, Air_CO2);
        Gb=BoundaryCound(1);
        Cb=BoundaryCound(2);
        Eb=BoundaryCound(3);
        
        % Compute stomatal conductance
        CalGs=ComputGsBallBerry(BallBerryIntercept,BallBerrySlope,WaterStressFactor,Eb,Cb,Gb,NetAssimilation,LeafTemperature,Air_CO2);
        Gs=CalGs(1);
        Ci=CalGs(2);
        % Compute energy balance
        CalLeafTemperature=ComputeEnergyBalance(Gb,Gs,NetAssimilation,LeafTemperature,WeatherRH,WeatherTemperature,Pressure,Boltzman,ConstantsCp,LatentHeatVaporization,Radiation_PAR,Radiation_NIR,Radiation_LW);
        LeafTemperature=CalLeafTemperature(1);
        Transpiration=CalLeafTemperature(2);
        % Compute convergence error
        ErrorLeafState_Ci = (PreviousLeafState_Ci - Ci) / Ci;
        ErrorLeafState_Gs = (PreviousLeafState_Gs - Gs) / Gs;
        ErrorLeafState_Temperature = (PreviousLeafState_Temperature - LeafTemperature) / LeafTemperature;

        %Check for oscillation and divergence to apply relaxation
        if (abs(Previous2Gs - Gs) < 0.01 && abs(ErrorLeafState_Gs) >= 0.001 && ErrorCount > 1) ||(abs(ErrorLeafState_Ci) > MaxError && ErrorCount > 1) || (abs(ErrorLeafState_Gs) > MaxError && ErrorCount > 1) ||(abs(ErrorLeafState_Temperature) > MaxError && ErrorCount > 1) % Divergence
            Relax = 0.5;
            % fprintf(LogOutputFile,"Relax");
        else if (Relax > 0.0)
            Relax = 0.0;
            end
        end
        %Apply relaxation due to oscillation and divergence
        Ci = Ci - Relax * (Ci - PreviousLeafState_Ci);
        if (Ci <= GammaStar)
            Ci = GammaStar; % ppm
        end
        Gs = Gs - Relax * (Gs - PreviousLeafState_Gs);
        if Gs < BallBerryIntercept
           Gs = BallBerryIntercept;
        end
        GrossAssimilation = GrossAssimilation - Relax *(GrossAssimilation- PreviousLeafMassFlux_GrossAssimilation);
        if (GrossAssimilation < 0.0)
            GrossAssimilation = 0.0;
            NetAssimilation = GrossAssimilation - Rd;
        end

        %Update values
        PreviousLeafState_Ci = Ci;
        Previous2Gs = PreviousLeafState_Gs; %Oscillation
        PreviousLeafState_Gs = Gs;
        PreviousLeafState_Temperature = LeafTemperature;
        PreviousLeafMassFlux_GrossAssimilation = GrossAssimilation;
        PreviousLeafMassFlux_NetAssimilation = NetAssimilation;
        ErrorCount = ErrorCount + 1;
    end
    
keySet={'Ci', 'NetAssimilation', 'GrossAssimilation', 'Gs', 'LeafTemperature', 'Transpiration'};
valueSet=[Ci, NetAssimilation, GrossAssimilation, Gs, LeafTemperature, Transpiration];
LeafA=containers.Map(keySet,valueSet);  

% LeafA(1)=Ci;
% LeafA(2)=NetAssimilation;
% LeafA(3)=Gs;
% LeafA(4)=LeafTemperature;
% LeafA(5)=Transpiration;
end
