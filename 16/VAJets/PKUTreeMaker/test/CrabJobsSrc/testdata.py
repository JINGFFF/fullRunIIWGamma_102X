infile =open('txtdata', 'r')
lines = infile.readlines()
#for i in range(0, 2):
for i in range(0, len(lines)):
    x = lines[i].split('%%%')
    name = x[0]
    pyfile = x[1]
    dataset = x[2]
    intxtfile = x[3]
    out_file_name = 'crab3_analysis' + name + '.py'
    outfile = open(out_file_name, 'w')
    outfile.write('from WMCore.Configuration import Configuration' + '\n')    
    outfile.write('config = Configuration()' + '\n')
    outfile.write('config.section_("General")' + '\n')
    outfile.write("config.General.requestName   = '" + name + "'" + '\n')
    outfile.write('config.General.transferLogs = True' + '\n')
    outfile.write('\n')

    outfile.write('config.section_("JobType")' + '\n')
    outfile.write("config.JobType.pluginName  = 'Analysis'" + '\n')
    outfile.write(intxtfile + '\n')
    outfile.write("config.JobType.psetName    = '" + pyfile +"'" + '\n')
    outfile.write('config.JobType.allowUndistributedCMSSW = True' + '\n')
    outfile.write('\n')

    outfile.write('config.section_("Data")' + '\n')
    outfile.write("config.Data.outputPrimaryDataset = 'VBS_WGAMMA_94X'" + '\n')
    outfile.write("config.Data.inputDataset = '" + dataset + "'"+ '\n')
    outfile.write("config.Data.inputDBS = 'global'" + '\n')
    outfile.write("config.Data.splitting = 'LumiBased'" + '\n')
    outfile.write('config.Data.unitsPerJob = 40' + '\n')
    outfile.write("config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'" + '\n')
    outfile.write('config.Data.publication = False' + '\n')
    outfile.write("config.Data.outputDatasetTag = '" + name + "'" + '\n')
    outfile.write('\n')

    outfile.write('config.section_("Site")' + '\n')
    outfile.write("config.Site.storageSite = 'T2_CN_Beijing'" + '\n')

    outfile.close()
    print(x[1])
infile.close()
