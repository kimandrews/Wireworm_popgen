#!/bin/bash

gatk --java-options "-Xmx900g" HaplotypeCaller  -R /mnt/local-scratch1/GRC_Projects/GRC_Wireworm/2019-04-09-WW.RAD_HiSeq_rev2/idx/idx.fasta -O ./05b-gatk/sub01_50reads_50map.vcf \
-I ./04-Mapped/CA02D_dedup.bam \
-I ./04-Mapped/CA13C_dedup.bam \
-I ./04-Mapped/CA13E_dedup.bam \
-I ./04-Mapped/CA14A_dedup.bam \
-I ./04-Mapped/CA14B_dedup.bam \
-I ./04-Mapped/CA14C_dedup.bam \
-I ./04-Mapped/CA14D_dedup.bam \
-I ./04-Mapped/CA17A_dedup.bam \
-I ./04-Mapped/CA17B_dedup.bam \
-I ./04-Mapped/CA17C_dedup.bam \
-I ./04-Mapped/CA18A_dedup.bam \
-I ./04-Mapped/CA18C_dedup.bam \
-I ./04-Mapped/CA18D_dedup.bam \
-I ./04-Mapped/CA18E_dedup.bam \
-I ./04-Mapped/ID001_dedup.bam \
-I ./04-Mapped/ID002_dedup.bam \
-I ./04-Mapped/ID003_dedup.bam \
-I ./04-Mapped/ID005_dedup.bam \
-I ./04-Mapped/ID006_dedup.bam \
-I ./04-Mapped/ID009_dedup.bam \
-I ./04-Mapped/ID010_dedup.bam \
-I ./04-Mapped/ID011_dedup.bam \
-I ./04-Mapped/ID014_dedup.bam \
-I ./04-Mapped/ID016_dedup.bam \
-I ./04-Mapped/ID017_dedup.bam \
-I ./04-Mapped/ID020_dedup.bam \
-I ./04-Mapped/ID021_dedup.bam \
-I ./04-Mapped/ID026_dedup.bam \
-I ./04-Mapped/ID029_dedup.bam \
-I ./04-Mapped/ID030_dedup.bam \
-I ./04-Mapped/ID031_dedup.bam \
-I ./04-Mapped/ID032_dedup.bam \
-I ./04-Mapped/ID033_dedup.bam \
-I ./04-Mapped/ID034_dedup.bam \
-I ./04-Mapped/ID035_dedup.bam \
-I ./04-Mapped/ID036_dedup.bam \
-I ./04-Mapped/ID037_dedup.bam \
-I ./04-Mapped/ID038_dedup.bam \
-I ./04-Mapped/ID039_dedup.bam \
-I ./04-Mapped/ID041_dedup.bam \
-I ./04-Mapped/ID042_dedup.bam \
-I ./04-Mapped/ID044_dedup.bam \
-I ./04-Mapped/ID046_dedup.bam \
-I ./04-Mapped/ID047_dedup.bam \
-I ./04-Mapped/ID054_dedup.bam \
-I ./04-Mapped/ID057_dedup.bam \
-I ./04-Mapped/ID061_dedup.bam \
-I ./04-Mapped/ID063_dedup.bam \
-I ./04-Mapped/ID064_dedup.bam \
-I ./04-Mapped/ID103_dedup.bam \
-I ./04-Mapped/ID105_dedup.bam \
-I ./04-Mapped/ID107_dedup.bam \
-I ./04-Mapped/ID109_dedup.bam \
-I ./04-Mapped/ID112_dedup.bam \
-I ./04-Mapped/ID115_dedup.bam \
-I ./04-Mapped/ID116_dedup.bam \
-I ./04-Mapped/ID118_dedup.bam \
-I ./04-Mapped/ID119_dedup.bam \
-I ./04-Mapped/ID123_dedup.bam \
-I ./04-Mapped/ID130_dedup.bam \
-I ./04-Mapped/ID131_dedup.bam \
-I ./04-Mapped/ID132_dedup.bam \
-I ./04-Mapped/ID133_dedup.bam \
-I ./04-Mapped/ID135_dedup.bam \
-I ./04-Mapped/LI001_dedup.bam \
-I ./04-Mapped/LI002_dedup.bam \
-I ./04-Mapped/LI005_dedup.bam \
-I ./04-Mapped/LI006_dedup.bam \
-I ./04-Mapped/LI016_dedup.bam \
-I ./04-Mapped/LI018_dedup.bam \
-I ./04-Mapped/LI028_dedup.bam \
-I ./04-Mapped/LI029_dedup.bam \
-I ./04-Mapped/LI030_dedup.bam \
-I ./04-Mapped/LI031_dedup.bam \
-I ./04-Mapped/LI032_dedup.bam \
-I ./04-Mapped/LI033_dedup.bam \
-I ./04-Mapped/LI034_dedup.bam \
-I ./04-Mapped/LI035_dedup.bam \
-I ./04-Mapped/LI036_dedup.bam \
-I ./04-Mapped/LI037_dedup.bam \
-I ./04-Mapped/LI038_dedup.bam \
-I ./04-Mapped/LI039_dedup.bam \
-I ./04-Mapped/LI041_dedup.bam \
-I ./04-Mapped/LI042_dedup.bam \
-I ./04-Mapped/LI043_dedup.bam \
-I ./04-Mapped/LI044_dedup.bam \
-I ./04-Mapped/LI047_dedup.bam \
-I ./04-Mapped/LI048_dedup.bam \
-I ./04-Mapped/LI051_dedup.bam \
-I ./04-Mapped/LI053_dedup.bam \
-I ./04-Mapped/LI054_dedup.bam \
-I ./04-Mapped/LI056_dedup.bam \
-I ./04-Mapped/LI057_dedup.bam \
-I ./04-Mapped/LI059_dedup.bam \
-I ./04-Mapped/LI060_dedup.bam \
-I ./04-Mapped/LI061_dedup.bam \
-I ./04-Mapped/LI062_dedup.bam \
-I ./04-Mapped/MT001_dedup.bam \
-I ./04-Mapped/MT004_dedup.bam \
-I ./04-Mapped/MT005_dedup.bam \
-I ./04-Mapped/MT007_dedup.bam \
-I ./04-Mapped/OR001_dedup.bam \
-I ./04-Mapped/Or001_dedup.bam \
-I ./04-Mapped/OR002_dedup.bam \
-I ./04-Mapped/Or002_dedup.bam \
-I ./04-Mapped/Or003_dedup.bam \
-I ./04-Mapped/OR003_dedup.bam \
-I ./04-Mapped/OR004_dedup.bam \
-I ./04-Mapped/Or004_dedup.bam \
-I ./04-Mapped/OR005_dedup.bam \
-I ./04-Mapped/Or005_dedup.bam \
-I ./04-Mapped/OR006_dedup.bam \
-I ./04-Mapped/Or006_dedup.bam \
-I ./04-Mapped/OR007_dedup.bam \
-I ./04-Mapped/Or007_dedup.bam \
-I ./04-Mapped/Or008_dedup.bam \
-I ./04-Mapped/OR008_dedup.bam \
-I ./04-Mapped/Or009_dedup.bam \
-I ./04-Mapped/OR009_dedup.bam \
-I ./04-Mapped/OR010_dedup.bam \
-I ./04-Mapped/Or010_dedup.bam \
-I ./04-Mapped/OR011_dedup.bam \
-I ./04-Mapped/Or011_dedup.bam \
-I ./04-Mapped/OR012_dedup.bam \
-I ./04-Mapped/Or012_dedup.bam \
-I ./04-Mapped/OR013_dedup.bam \
-I ./04-Mapped/Or013_dedup.bam \
-I ./04-Mapped/Or014_dedup.bam \
-I ./04-Mapped/OR014_dedup.bam \
-I ./04-Mapped/Or015_dedup.bam \
-I ./04-Mapped/OR015_dedup.bam \
-I ./04-Mapped/OR016_dedup.bam \
-I ./04-Mapped/Or016_dedup.bam \
-I ./04-Mapped/OR017_dedup.bam \
-I ./04-Mapped/Or017_dedup.bam \
-I ./04-Mapped/OR018_dedup.bam \
-I ./04-Mapped/OR019_dedup.bam \
-I ./04-Mapped/OR020_dedup.bam \
-I ./04-Mapped/OR021_dedup.bam \
-I ./04-Mapped/OR022_dedup.bam \
-I ./04-Mapped/OR023_dedup.bam \
-I ./04-Mapped/OR024_dedup.bam \
-I ./04-Mapped/OR025_dedup.bam \
-I ./04-Mapped/OR026_dedup.bam \
-I ./04-Mapped/OR027_dedup.bam \
-I ./04-Mapped/OR028_dedup.bam \
-I ./04-Mapped/OR029_dedup.bam \
-I ./04-Mapped/OR031_dedup.bam \
-I ./04-Mapped/OR033_dedup.bam \
-I ./04-Mapped/OR034_dedup.bam \
-I ./04-Mapped/OR036_dedup.bam \
-I ./04-Mapped/OR037_dedup.bam \
-I ./04-Mapped/OR038_dedup.bam \
-I ./04-Mapped/OR039_dedup.bam \
-I ./04-Mapped/OR040_dedup.bam \
-I ./04-Mapped/OR041_dedup.bam \
-I ./04-Mapped/OR043_dedup.bam \
-I ./04-Mapped/OR044_dedup.bam \
-I ./04-Mapped/OR045_dedup.bam \
-I ./04-Mapped/OR046_dedup.bam \
-I ./04-Mapped/OR047_dedup.bam \
-I ./04-Mapped/OR049_dedup.bam \
-I ./04-Mapped/OR050_dedup.bam \
-I ./04-Mapped/OR051_dedup.bam \
-I ./04-Mapped/OR052_dedup.bam \
-I ./04-Mapped/OR053_dedup.bam \
-I ./04-Mapped/OR054_dedup.bam \
-I ./04-Mapped/OR055_dedup.bam \
-I ./04-Mapped/OR056_dedup.bam \
-I ./04-Mapped/OR057_dedup.bam \
-I ./04-Mapped/OR058_dedup.bam \
-I ./04-Mapped/OR059_dedup.bam \
-I ./04-Mapped/OR060_dedup.bam \
-I ./04-Mapped/OR062_dedup.bam \
-I ./04-Mapped/OR063_dedup.bam \
-I ./04-Mapped/OR064_dedup.bam \
-I ./04-Mapped/OR065_dedup.bam \
-I ./04-Mapped/OR066_dedup.bam \
-I ./04-Mapped/OR067_dedup.bam \
-I ./04-Mapped/OR068_dedup.bam \
-I ./04-Mapped/OR069_dedup.bam \
-I ./04-Mapped/OR070_dedup.bam \
-I ./04-Mapped/OR071_dedup.bam \
-I ./04-Mapped/OR072_dedup.bam \
-I ./04-Mapped/OR073_dedup.bam \
-I ./04-Mapped/OR074_dedup.bam \
-I ./04-Mapped/OR075_dedup.bam \
-I ./04-Mapped/OR076_dedup.bam \
-I ./04-Mapped/OR077_dedup.bam \
-I ./04-Mapped/OR078_dedup.bam \
-I ./04-Mapped/OR079_dedup.bam \
-I ./04-Mapped/OR080_dedup.bam \
-I ./04-Mapped/OR081_dedup.bam \
-I ./04-Mapped/OR082_dedup.bam \
-I ./04-Mapped/OR083_dedup.bam \
-I ./04-Mapped/OR084_dedup.bam \
-I ./04-Mapped/OR085_dedup.bam \
-I ./04-Mapped/OR086_dedup.bam \
-I ./04-Mapped/OR087_dedup.bam \
-I ./04-Mapped/OR088_dedup.bam \
-I ./04-Mapped/OR089_dedup.bam \
-I ./04-Mapped/OR090_dedup.bam \
-I ./04-Mapped/OR091_dedup.bam \
-I ./04-Mapped/OR092_dedup.bam \
-I ./04-Mapped/OR093_dedup.bam \
-I ./04-Mapped/OR094_dedup.bam \
-I ./04-Mapped/OR095_dedup.bam \
-I ./04-Mapped/OR096_dedup.bam \
-I ./04-Mapped/OR097_dedup.bam \
-I ./04-Mapped/OR098_dedup.bam \
-I ./04-Mapped/OR099_dedup.bam \
-I ./04-Mapped/OR100_dedup.bam \
-I ./04-Mapped/OR101_dedup.bam \
-I ./04-Mapped/OR102_dedup.bam \
-I ./04-Mapped/OR103_dedup.bam \
-I ./04-Mapped/OR104_dedup.bam \
-I ./04-Mapped/OR105_dedup.bam \
-I ./04-Mapped/OR106_dedup.bam \
-I ./04-Mapped/OR107_dedup.bam \
-I ./04-Mapped/OR108_dedup.bam \
-I ./04-Mapped/OR109_dedup.bam \
-I ./04-Mapped/OR110_dedup.bam \
-I ./04-Mapped/OR111_dedup.bam \
-I ./04-Mapped/OR112_dedup.bam \
-I ./04-Mapped/OR113_dedup.bam \
-I ./04-Mapped/OR114_dedup.bam \
-I ./04-Mapped/OR115_dedup.bam \
-I ./04-Mapped/OR116_dedup.bam \
-I ./04-Mapped/OR117_dedup.bam \
-I ./04-Mapped/OR118_dedup.bam \
-I ./04-Mapped/OR119_dedup.bam \
-I ./04-Mapped/OR120_dedup.bam \
-I ./04-Mapped/OR122_dedup.bam \
-I ./04-Mapped/OR123_dedup.bam \
-I ./04-Mapped/OR124_dedup.bam \
-I ./04-Mapped/OR125_dedup.bam \
-I ./04-Mapped/OR127_dedup.bam \
-I ./04-Mapped/OR128_dedup.bam \
-I ./04-Mapped/OR129_dedup.bam \
-I ./04-Mapped/OR130_dedup.bam \
-I ./04-Mapped/OR131_dedup.bam \
-I ./04-Mapped/OR132_dedup.bam \
-I ./04-Mapped/OR134_dedup.bam \
-I ./04-Mapped/OR135_dedup.bam \
-I ./04-Mapped/OR136_dedup.bam \
-I ./04-Mapped/OR137_dedup.bam \
-I ./04-Mapped/OR138_dedup.bam \
-I ./04-Mapped/OR140_dedup.bam \
-I ./04-Mapped/OR141_dedup.bam \
-I ./04-Mapped/OR142_dedup.bam \
-I ./04-Mapped/OR143_dedup.bam \
-I ./04-Mapped/OR144_dedup.bam \
-I ./04-Mapped/OR145_dedup.bam \
-I ./04-Mapped/OR146_dedup.bam \
-I ./04-Mapped/OR147_dedup.bam \
-I ./04-Mapped/OR148_dedup.bam \
-I ./04-Mapped/OR149_dedup.bam \
-I ./04-Mapped/OR150_dedup.bam \
-I ./04-Mapped/OR151_dedup.bam \
-I ./04-Mapped/OR153_dedup.bam \
-I ./04-Mapped/OR154_dedup.bam \
-I ./04-Mapped/OR155_dedup.bam \
-I ./04-Mapped/OR156_dedup.bam \
-I ./04-Mapped/OR157_dedup.bam \
-I ./04-Mapped/OR158_dedup.bam \
-I ./04-Mapped/OR159_dedup.bam \
-I ./04-Mapped/OR160_dedup.bam \
-I ./04-Mapped/OR161_dedup.bam \
-I ./04-Mapped/OR162_dedup.bam \
-I ./04-Mapped/OR163_dedup.bam \
-I ./04-Mapped/OR164_dedup.bam \
-I ./04-Mapped/OR165_dedup.bam \
-I ./04-Mapped/OR166_dedup.bam \
-I ./04-Mapped/OR167_dedup.bam \
-I ./04-Mapped/OR168_dedup.bam \
-I ./04-Mapped/OR169_dedup.bam \
-I ./04-Mapped/OR170_dedup.bam \
-I ./04-Mapped/OR171_dedup.bam \
-I ./04-Mapped/OR172_dedup.bam \
-I ./04-Mapped/OR173_dedup.bam \
-I ./04-Mapped/OR174_dedup.bam \
-I ./04-Mapped/OR175_dedup.bam \
-I ./04-Mapped/OR176_dedup.bam \
-I ./04-Mapped/OR177_dedup.bam \
-I ./04-Mapped/WA001_dedup.bam \
-I ./04-Mapped/WA002_dedup.bam \
-I ./04-Mapped/WA003_dedup.bam \
-I ./04-Mapped/WA004_dedup.bam \
-I ./04-Mapped/WA005_dedup.bam \
-I ./04-Mapped/WA006_dedup.bam \
-I ./04-Mapped/WA007_dedup.bam \
-I ./04-Mapped/WA008_dedup.bam \
-I ./04-Mapped/WA009_dedup.bam \
-I ./04-Mapped/WA010_dedup.bam \
-I ./04-Mapped/WA011_dedup.bam \
-I ./04-Mapped/WA012_dedup.bam \
-I ./04-Mapped/WA013_dedup.bam \
-I ./04-Mapped/WA014_dedup.bam \
-I ./04-Mapped/WA015_dedup.bam \
-I ./04-Mapped/WA016_dedup.bam \
-I ./04-Mapped/WA017_dedup.bam \
-I ./04-Mapped/WA018_dedup.bam \
-I ./04-Mapped/WA019_dedup.bam \
-I ./04-Mapped/WA020_dedup.bam \
-I ./04-Mapped/WA021_dedup.bam \
-I ./04-Mapped/WA022_dedup.bam \
-I ./04-Mapped/WA023_dedup.bam \
-I ./04-Mapped/WA024_dedup.bam \
-I ./04-Mapped/WA025_dedup.bam \
-I ./04-Mapped/WA026_dedup.bam \
-I ./04-Mapped/WA028_dedup.bam \
-I ./04-Mapped/WA030_dedup.bam \
-I ./04-Mapped/WA031_dedup.bam \
-I ./04-Mapped/WA032_dedup.bam \
-I ./04-Mapped/WA033_dedup.bam \
-I ./04-Mapped/WA034_dedup.bam \
-I ./04-Mapped/WA048_dedup.bam \
-I ./04-Mapped/WA049_dedup.bam \
-I ./04-Mapped/WA050_dedup.bam \
-I ./04-Mapped/WA051_dedup.bam \
-I ./04-Mapped/WA053_dedup.bam \
-I ./04-Mapped/WA063_dedup.bam \
-I ./04-Mapped/WA065_dedup.bam