function LeafAssim=LeafWrap2(Env, Einput, Eio)
  global pcfactor;  
  pcfactor=1.0/Env('ProteinTotalRatio');
  % Rearrange Einput to use the same order as sample
  inE=importdata('MeM_input5_0.txt');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Input conversion
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
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Einput=inE.data;
  Einput=ExpValue;
  Edata=importdata('Einput7.txt');
  Eio=Edata.data(:,1);
  MetaOnly=0;% if MetaOnly=1 run Metabolic model
  
  WeatherTemp=25;
  Air_CO2=400;
  WeatherRH=0.6;
  WeatherWind=5;
  Radiation_PAR=470*0.85*0.85;%%%%
  Convert=1E6/(2.35E5); %Convert W m^{-2} to u moles m^{-2} s^{-1}
  Radiation_NIR=0;
  Radiation_LW=0;
  PhotosynthesisType=1.1;
  Vcmax25=100;
  Jmax25=200;
  GRNC=1;

  if MetaOnly==1
    CO2i=Air_CO2*0.7; % intercellular CO2 
    PPFDi=Radiation_PAR*Convert;
    NetAssimilation=EPS_Drive_GRNs(Einput,CO2i,PPFDi,WeatherTemp,GRNC,0,Eio);
  else
    LeafResult=Leaf(Env('WeatherRH'),Env('WeatherTemperature'),Env('Air_CO2'),...
        Env('WeatherWind'),Env('Radiation_PAR'),Env('Radiation_NIR'),...
        Env('Radiation_LW'),Env('PhotosynthesisType'),Env('Vcmax25'),...
        Env('Jmax25'),Env('GRNC'), ExpValue, Eio{:, 2});
    Ci=LeafResult(1);
    NetAssimilation=LeafResult(2);
    Gs=LeafResult(3);
    LeafTemperature=LeafResult(4);
    Transpirationi=LeafResult(5);
  end
  
  LeafAssim = NetAssimilation;
  
end