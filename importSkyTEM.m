function YellowstoneNP2016NorrisProcessedDataModels = importSkyTEM(filename, dataLines)
%IMPORTFILE Import data from a text file
%  YELLOWSTONENP2016NORRISPROCESSEDDATAMODELS = IMPORTFILE(FILENAME)
%  reads data from text file FILENAME for the default selection.
%  Returns the numeric data.
%
%  YELLOWSTONENP2016NORRISPROCESSEDDATAMODELS = IMPORTFILE(FILE,
%  DATALINES) reads data for the specified row interval(s) of text file
%  FILENAME. Specify DATALINES as a positive scalar integer or a N-by-2
%  array of positive scalar integers for dis-contiguous row intervals.
%
%  Example:
%  YellowstoneNP2016NorrisProcessedDataModels = importfile("C:\Users\aparseki\Dropbox\YNP\USGS_YNP_SkyTEM\AirborneElectro\YellowstoneNP_2016_Norris_ProcessedData_Models.csv", [1, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 30-Nov-2023 14:16:16

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [1, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 442);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Flight", "Var3", "E_UTM12N", "N_UTM12N", "Ytopbathv2", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "RHO_I0", "RHO_I1", "RHO_I2", "RHO_I3", "RHO_I4", "RHO_I5", "RHO_I6", "RHO_I7", "RHO_I8", "RHO_I9", "RHO_I10", "RHO_I11", "RHO_I12", "RHO_I13", "RHO_I14", "RHO_I15", "RHO_I16", "RHO_I17", "RHO_I18", "RHO_I19", "RHO_I20", "RHO_I21", "RHO_I22", "RHO_I23", "RHO_I24", "RHO_I25", "RHO_I26", "RHO_I27", "RHO_I28", "RHO_I29", "RHO_I_STD0", "RHO_I_STD1", "RHO_I_STD2", "RHO_I_STD3", "RHO_I_STD4", "RHO_I_STD5", "RHO_I_STD6", "RHO_I_STD7", "RHO_I_STD8", "RHO_I_STD9", "RHO_I_STD10", "RHO_I_STD11", "RHO_I_STD12", "RHO_I_STD13", "RHO_I_STD14", "RHO_I_STD15", "RHO_I_STD16", "RHO_I_STD17", "RHO_I_STD18", "RHO_I_STD19", "RHO_I_STD20", "RHO_I_STD21", "RHO_I_STD22", "RHO_I_STD23", "RHO_I_STD24", "RHO_I_STD25", "RHO_I_STD26", "RHO_I_STD27", "RHO_I_STD28", "RHO_I_STD29", "DEP_TOP0", "DEP_TOP1", "DEP_TOP2", "DEP_TOP3", "DEP_TOP4", "DEP_TOP5", "DEP_TOP6", "DEP_TOP7", "DEP_TOP8", "DEP_TOP9", "DEP_TOP10", "DEP_TOP11", "DEP_TOP12", "DEP_TOP13", "DEP_TOP14", "DEP_TOP15", "DEP_TOP16", "DEP_TOP17", "DEP_TOP18", "DEP_TOP19", "DEP_TOP20", "DEP_TOP21", "DEP_TOP22", "DEP_TOP23", "DEP_TOP24", "DEP_TOP25", "DEP_TOP26", "DEP_TOP27", "DEP_TOP28", "DEP_TOP29", "DEP_BOT0", "DEP_BOT1", "DEP_BOT2", "DEP_BOT3", "DEP_BOT4", "DEP_BOT5", "DEP_BOT6", "DEP_BOT7", "DEP_BOT8", "DEP_BOT9", "DEP_BOT10", "DEP_BOT11", "DEP_BOT12", "DEP_BOT13", "DEP_BOT14", "DEP_BOT15", "DEP_BOT16", "DEP_BOT17", "DEP_BOT18", "DEP_BOT19", "DEP_BOT20", "DEP_BOT21", "DEP_BOT22", "DEP_BOT23", "DEP_BOT24", "DEP_BOT25", "DEP_BOT26", "DEP_BOT27", "DEP_BOT28", "DEP_BOT29", "THK0", "THK1", "THK2", "THK3", "THK4", "THK5", "THK6", "THK7", "THK8", "THK9", "THK10", "THK11", "THK12", "THK13", "THK14", "THK15", "THK16", "THK17", "THK18", "THK19", "THK20", "THK21", "THK22", "THK23", "THK24", "THK25", "THK26", "THK27", "THK28", "THK_STD0", "THK_STD1", "THK_STD2", "THK_STD3", "THK_STD4", "THK_STD5", "THK_STD6", "THK_STD7", "THK_STD8", "THK_STD9", "THK_STD10", "THK_STD11", "THK_STD12", "THK_STD13", "THK_STD14", "THK_STD15", "THK_STD16", "THK_STD17", "THK_STD18", "THK_STD19", "THK_STD20", "THK_STD21", "THK_STD22", "THK_STD23", "THK_STD24", "THK_STD25", "THK_STD26", "THK_STD27", "THK_STD28", "DEP_BOT_STD0", "DEP_BOT_STD1", "DEP_BOT_STD2", "DEP_BOT_STD3", "DEP_BOT_STD4", "DEP_BOT_STD5", "DEP_BOT_STD6", "DEP_BOT_STD7", "DEP_BOT_STD8", "DEP_BOT_STD9", "DEP_BOT_STD10", "DEP_BOT_STD11", "DEP_BOT_STD12", "DEP_BOT_STD13", "DEP_BOT_STD14", "DEP_BOT_STD15", "DEP_BOT_STD16", "DEP_BOT_STD17", "DEP_BOT_STD18", "DEP_BOT_STD19", "DEP_BOT_STD20", "DEP_BOT_STD21", "DEP_BOT_STD22", "DEP_BOT_STD23", "DEP_BOT_STD24", "DEP_BOT_STD25", "DEP_BOT_STD26", "DEP_BOT_STD27", "DEP_BOT_STD28", "DOI_STANDARD", "Var223", "Var224", "Var225", "Var226", "Var227", "Var228", "Var229", "Var230", "Var231", "Var232", "Var233", "Var234", "Var235", "Var236", "Var237", "Var238", "Var239", "Var240", "Var241", "Var242", "Var243", "Var244", "Var245", "Var246", "Var247", "Var248", "Var249", "Var250", "Var251", "Var252", "Var253", "Var254", "Var255", "Var256", "Var257", "Var258", "Var259", "Var260", "Var261", "Var262", "Var263", "Var264", "Var265", "Var266", "Var267", "Var268", "Var269", "Var270", "Var271", "Var272", "Var273", "Var274", "Var275", "Var276", "Var277", "Var278", "Var279", "Var280", "Var281", "Var282", "Var283", "Var284", "Var285", "Var286", "Var287", "Var288", "Var289", "Var290", "Var291", "Var292", "Var293", "Var294", "Var295", "Var296", "Var297", "Var298", "Var299", "Var300", "Var301", "Var302", "Var303", "Var304", "Var305", "Var306", "Var307", "Var308", "Var309", "Var310", "Var311", "Var312", "Var313", "Var314", "Var315", "Var316", "Var317", "Var318", "Var319", "Var320", "Var321", "Var322", "Var323", "Var324", "Var325", "Var326", "Var327", "Var328", "Var329", "Var330", "Var331", "Var332", "Var333", "Var334", "Var335", "Var336", "Var337", "Var338", "Var339", "Var340", "Var341", "Var342", "Var343", "Var344", "Var345", "Var346", "Var347", "Var348", "Var349", "Var350", "Var351", "Var352", "Var353", "Var354", "Var355", "Var356", "Var357", "Var358", "Var359", "Var360", "Var361", "Var362", "Var363", "Var364", "Var365", "Var366", "Var367", "Var368", "Var369", "Var370", "Var371", "Var372", "Var373", "Var374", "Var375", "Var376", "Var377", "Var378", "Var379", "Var380", "Var381", "Var382", "Var383", "Var384", "Var385", "Var386", "Var387", "Var388", "Var389", "Var390", "Var391", "Var392", "Var393", "Var394", "Var395", "Var396", "Var397", "Var398", "Var399", "Var400", "Var401", "Var402", "Var403", "Var404", "Var405", "Var406", "Var407", "Var408", "Var409", "Var410", "Var411", "Var412", "Var413", "Var414", "Var415", "Var416", "Var417", "Var418", "Var419", "Var420", "Var421", "Var422", "Var423", "Var424", "Var425", "Var426", "Var427", "Var428", "Var429", "Var430", "Var431", "Var432", "Var433", "Var434", "Var435", "Var436", "Var437", "Var438", "Var439", "Var440", "Var441", "Var442"];
opts.SelectedVariableNames = ["Flight", "E_UTM12N", "N_UTM12N", "Ytopbathv2", "RHO_I0", "RHO_I1", "RHO_I2", "RHO_I3", "RHO_I4", "RHO_I5", "RHO_I6", "RHO_I7", "RHO_I8", "RHO_I9", "RHO_I10", "RHO_I11", "RHO_I12", "RHO_I13", "RHO_I14", "RHO_I15", "RHO_I16", "RHO_I17", "RHO_I18", "RHO_I19", "RHO_I20", "RHO_I21", "RHO_I22", "RHO_I23", "RHO_I24", "RHO_I25", "RHO_I26", "RHO_I27", "RHO_I28", "RHO_I29", "RHO_I_STD0", "RHO_I_STD1", "RHO_I_STD2", "RHO_I_STD3", "RHO_I_STD4", "RHO_I_STD5", "RHO_I_STD6", "RHO_I_STD7", "RHO_I_STD8", "RHO_I_STD9", "RHO_I_STD10", "RHO_I_STD11", "RHO_I_STD12", "RHO_I_STD13", "RHO_I_STD14", "RHO_I_STD15", "RHO_I_STD16", "RHO_I_STD17", "RHO_I_STD18", "RHO_I_STD19", "RHO_I_STD20", "RHO_I_STD21", "RHO_I_STD22", "RHO_I_STD23", "RHO_I_STD24", "RHO_I_STD25", "RHO_I_STD26", "RHO_I_STD27", "RHO_I_STD28", "RHO_I_STD29", "DEP_TOP0", "DEP_TOP1", "DEP_TOP2", "DEP_TOP3", "DEP_TOP4", "DEP_TOP5", "DEP_TOP6", "DEP_TOP7", "DEP_TOP8", "DEP_TOP9", "DEP_TOP10", "DEP_TOP11", "DEP_TOP12", "DEP_TOP13", "DEP_TOP14", "DEP_TOP15", "DEP_TOP16", "DEP_TOP17", "DEP_TOP18", "DEP_TOP19", "DEP_TOP20", "DEP_TOP21", "DEP_TOP22", "DEP_TOP23", "DEP_TOP24", "DEP_TOP25", "DEP_TOP26", "DEP_TOP27", "DEP_TOP28", "DEP_TOP29", "DEP_BOT0", "DEP_BOT1", "DEP_BOT2", "DEP_BOT3", "DEP_BOT4", "DEP_BOT5", "DEP_BOT6", "DEP_BOT7", "DEP_BOT8", "DEP_BOT9", "DEP_BOT10", "DEP_BOT11", "DEP_BOT12", "DEP_BOT13", "DEP_BOT14", "DEP_BOT15", "DEP_BOT16", "DEP_BOT17", "DEP_BOT18", "DEP_BOT19", "DEP_BOT20", "DEP_BOT21", "DEP_BOT22", "DEP_BOT23", "DEP_BOT24", "DEP_BOT25", "DEP_BOT26", "DEP_BOT27", "DEP_BOT28", "DEP_BOT29", "THK0", "THK1", "THK2", "THK3", "THK4", "THK5", "THK6", "THK7", "THK8", "THK9", "THK10", "THK11", "THK12", "THK13", "THK14", "THK15", "THK16", "THK17", "THK18", "THK19", "THK20", "THK21", "THK22", "THK23", "THK24", "THK25", "THK26", "THK27", "THK28", "THK_STD0", "THK_STD1", "THK_STD2", "THK_STD3", "THK_STD4", "THK_STD5", "THK_STD6", "THK_STD7", "THK_STD8", "THK_STD9", "THK_STD10", "THK_STD11", "THK_STD12", "THK_STD13", "THK_STD14", "THK_STD15", "THK_STD16", "THK_STD17", "THK_STD18", "THK_STD19", "THK_STD20", "THK_STD21", "THK_STD22", "THK_STD23", "THK_STD24", "THK_STD25", "THK_STD26", "THK_STD27", "THK_STD28", "DEP_BOT_STD0", "DEP_BOT_STD1", "DEP_BOT_STD2", "DEP_BOT_STD3", "DEP_BOT_STD4", "DEP_BOT_STD5", "DEP_BOT_STD6", "DEP_BOT_STD7", "DEP_BOT_STD8", "DEP_BOT_STD9", "DEP_BOT_STD10", "DEP_BOT_STD11", "DEP_BOT_STD12", "DEP_BOT_STD13", "DEP_BOT_STD14", "DEP_BOT_STD15", "DEP_BOT_STD16", "DEP_BOT_STD17", "DEP_BOT_STD18", "DEP_BOT_STD19", "DEP_BOT_STD20", "DEP_BOT_STD21", "DEP_BOT_STD22", "DEP_BOT_STD23", "DEP_BOT_STD24", "DEP_BOT_STD25", "DEP_BOT_STD26", "DEP_BOT_STD27", "DEP_BOT_STD28", "DOI_STANDARD"];
opts.VariableTypes = ["string", "double", "string", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var3", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var223", "Var224", "Var225", "Var226", "Var227", "Var228", "Var229", "Var230", "Var231", "Var232", "Var233", "Var234", "Var235", "Var236", "Var237", "Var238", "Var239", "Var240", "Var241", "Var242", "Var243", "Var244", "Var245", "Var246", "Var247", "Var248", "Var249", "Var250", "Var251", "Var252", "Var253", "Var254", "Var255", "Var256", "Var257", "Var258", "Var259", "Var260", "Var261", "Var262", "Var263", "Var264", "Var265", "Var266", "Var267", "Var268", "Var269", "Var270", "Var271", "Var272", "Var273", "Var274", "Var275", "Var276", "Var277", "Var278", "Var279", "Var280", "Var281", "Var282", "Var283", "Var284", "Var285", "Var286", "Var287", "Var288", "Var289", "Var290", "Var291", "Var292", "Var293", "Var294", "Var295", "Var296", "Var297", "Var298", "Var299", "Var300", "Var301", "Var302", "Var303", "Var304", "Var305", "Var306", "Var307", "Var308", "Var309", "Var310", "Var311", "Var312", "Var313", "Var314", "Var315", "Var316", "Var317", "Var318", "Var319", "Var320", "Var321", "Var322", "Var323", "Var324", "Var325", "Var326", "Var327", "Var328", "Var329", "Var330", "Var331", "Var332", "Var333", "Var334", "Var335", "Var336", "Var337", "Var338", "Var339", "Var340", "Var341", "Var342", "Var343", "Var344", "Var345", "Var346", "Var347", "Var348", "Var349", "Var350", "Var351", "Var352", "Var353", "Var354", "Var355", "Var356", "Var357", "Var358", "Var359", "Var360", "Var361", "Var362", "Var363", "Var364", "Var365", "Var366", "Var367", "Var368", "Var369", "Var370", "Var371", "Var372", "Var373", "Var374", "Var375", "Var376", "Var377", "Var378", "Var379", "Var380", "Var381", "Var382", "Var383", "Var384", "Var385", "Var386", "Var387", "Var388", "Var389", "Var390", "Var391", "Var392", "Var393", "Var394", "Var395", "Var396", "Var397", "Var398", "Var399", "Var400", "Var401", "Var402", "Var403", "Var404", "Var405", "Var406", "Var407", "Var408", "Var409", "Var410", "Var411", "Var412", "Var413", "Var414", "Var415", "Var416", "Var417", "Var418", "Var419", "Var420", "Var421", "Var422", "Var423", "Var424", "Var425", "Var426", "Var427", "Var428", "Var429", "Var430", "Var431", "Var432", "Var433", "Var434", "Var435", "Var436", "Var437", "Var438", "Var439", "Var440", "Var441", "Var442"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var3", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var223", "Var224", "Var225", "Var226", "Var227", "Var228", "Var229", "Var230", "Var231", "Var232", "Var233", "Var234", "Var235", "Var236", "Var237", "Var238", "Var239", "Var240", "Var241", "Var242", "Var243", "Var244", "Var245", "Var246", "Var247", "Var248", "Var249", "Var250", "Var251", "Var252", "Var253", "Var254", "Var255", "Var256", "Var257", "Var258", "Var259", "Var260", "Var261", "Var262", "Var263", "Var264", "Var265", "Var266", "Var267", "Var268", "Var269", "Var270", "Var271", "Var272", "Var273", "Var274", "Var275", "Var276", "Var277", "Var278", "Var279", "Var280", "Var281", "Var282", "Var283", "Var284", "Var285", "Var286", "Var287", "Var288", "Var289", "Var290", "Var291", "Var292", "Var293", "Var294", "Var295", "Var296", "Var297", "Var298", "Var299", "Var300", "Var301", "Var302", "Var303", "Var304", "Var305", "Var306", "Var307", "Var308", "Var309", "Var310", "Var311", "Var312", "Var313", "Var314", "Var315", "Var316", "Var317", "Var318", "Var319", "Var320", "Var321", "Var322", "Var323", "Var324", "Var325", "Var326", "Var327", "Var328", "Var329", "Var330", "Var331", "Var332", "Var333", "Var334", "Var335", "Var336", "Var337", "Var338", "Var339", "Var340", "Var341", "Var342", "Var343", "Var344", "Var345", "Var346", "Var347", "Var348", "Var349", "Var350", "Var351", "Var352", "Var353", "Var354", "Var355", "Var356", "Var357", "Var358", "Var359", "Var360", "Var361", "Var362", "Var363", "Var364", "Var365", "Var366", "Var367", "Var368", "Var369", "Var370", "Var371", "Var372", "Var373", "Var374", "Var375", "Var376", "Var377", "Var378", "Var379", "Var380", "Var381", "Var382", "Var383", "Var384", "Var385", "Var386", "Var387", "Var388", "Var389", "Var390", "Var391", "Var392", "Var393", "Var394", "Var395", "Var396", "Var397", "Var398", "Var399", "Var400", "Var401", "Var402", "Var403", "Var404", "Var405", "Var406", "Var407", "Var408", "Var409", "Var410", "Var411", "Var412", "Var413", "Var414", "Var415", "Var416", "Var417", "Var418", "Var419", "Var420", "Var421", "Var422", "Var423", "Var424", "Var425", "Var426", "Var427", "Var428", "Var429", "Var430", "Var431", "Var432", "Var433", "Var434", "Var435", "Var436", "Var437", "Var438", "Var439", "Var440", "Var441", "Var442"], "EmptyFieldRule", "auto");

% Import the data
YellowstoneNP2016NorrisProcessedDataModels = readtable(filename, opts);

%% Convert to output type
YellowstoneNP2016NorrisProcessedDataModels = table2array(YellowstoneNP2016NorrisProcessedDataModels);
end