library(readr)
library(plyr)
library(tidyverse)
library(data.table)
library(utils)
library(tcltk)
library(ggrepel)
library(ggsci)
library(strex)
library(varhandle)

###

#compoundID looks in the first line of the .txt file - make sure file name is "*_MSMSINFUSION_compoundname_*", and is mirrored in the first line of the .txt file
#compoundID = str_split_fixed(str_split_fixed(as.character(all_files[[1]][[1]][1]), "MSMSINFUSION_*", 2)[[2]], "_", 2)[[1]]
#alternatively, just use readline()
compoundID = readline(prompt="Enter compound ID: ")

rm(list=setdiff(ls(),"compoundID"))

#import files (shift+click, ctrl+click)

files = tk_choose.files()
all_files = lapply(files, function(x){read_tsv(file = x)})

#or, by directory:
#files = dir_ls(path = choose.dir(), glob = "*txt")
#read_tsv(files, col_names = F)
#all_files = lapply(files, function(x){read.table(file = x, sep = '\t', header = F)})

instrument = c("Unknown instrument","Unknown","Waters","Waters instrument", "LTQ Orbitrap Discovery - Ion Trap", "ITMS", "LTQ Orbitrap Discovery - FTMS", "FTMS")
model = 1
c_mass = 166.08 # this will likely need to be readline() in future
c_charge = c(1,0.9,0.85,0.8,0.75,0.75) #charge adjustment factors for thermo LTQ HCD collision energy calculation = +1 = 1, +2 = 0.9, +3 = 0.85 etc
charge = 1

if (any(str_detect(all_files[[1]],"cvParam: Waters"))) {model = 3} else 
  if (any(str_detect(all_files[[1]],"cvParam: LTQ Orbitrap Discovery"))) 
    {
    if (any(str_detect(all_files[[1]],"cvParam: filter string, ITMS"))) {model = 5} else
    if (any(str_detect(all_files[[1]],"cvParam: filter string, FTMS"))) {model = 7}
    } else model = 1

###

CElist = list("mass")
exptlist = list()
for (j in 1:length(all_files))
{
  my_data = all_files[[j]]
  colnames(my_data) = "column"
  #Get scan width and collision energy for expt
  combinewidth = as.numeric(str_split_fixed(str_split_fixed(as.character(my_data[str_detect(my_data[[1]],"spectrumList "),]), "\\(", 2)[[2]], " ", 2)[[1]])
  CE = round(as.numeric(str_split_fixed(str_split_fixed(as.character(my_data[str_detect(my_data[[1]],"cvParam: collision energy, "),]), "\\, electronvolt", 2)[[1]], "\\, ", 2)[[2]]), digits = 1)
  if (model == 7) {CE = round_any(c_charge[charge]*c_mass/500*CE,0.1)}


  intensities = my_data %>% 
    filter(str_detect(column,"binary: ")) %>%
    slice(1:(n()-3))

  scanlist = list()

  for (i in 1:combinewidth)
  {
    temp = data.frame(round_any(str_extract_numbers(slice(intensities,(2*i-1):(2*i))[[1]][1],decimals = T, sci = T)[[1]][-1], ifelse(model == 5,1,0.05), f = round), str_extract_numbers(slice(intensities,(2*i-1):(2*i))[[1]][2],decimals = T, sci = T)[[1]][-1])
    colnames(temp) = c("Mass","Intensity")
    temp = temp %>%
     #filter(temp[2] > 0) %>%
     group_by_at(1) %>%
     summarise(across(everything(),if (model == 7 | model == 5){mean}else{
       sum
       }
       ))
   
   scanlist = append(scanlist, list(temp))
  }

  CElist = append(CElist, as.character(CE))
  exptlist = append(exptlist, list(summarise(group_by_at(rbindlist(scanlist),"Mass"), across(everything(),sum))))
  significant = ifelse(exptlist[[j]][["Intensity"]]>max(exptlist[[j]][["Intensity"]],na.rm=T)/20, paste(as.character(exptlist[[j]][["Mass"]]), as.character(round_any(exptlist[[j]][["Intensity"]],1)), sep = " ; "),"")
  spectrum = ggplot(data=exptlist[[j]],mapping = aes(x=Mass,y = Intensity, ymax=Intensity ,ymin=0)) + 
    geom_linerange() +
    xlim(50,170) + ylim(0,(max(exptlist[[j]][["Intensity"]],na.rm=T)+1000)) + 
    geom_text_repel(point.size = NA, min.segment.length=0, aes(label=significant, vjust = -0.5)) +
    labs(x = "m/z", 
         y = "Intensity", 
         title = str_wrap(paste(instrument[model], "MS/MS spectrum for ", compoundID, " @ ", CE), width = 60))
  
  png(filename=paste(instrument[model+1], compoundID, CE, "MSMS.png", sep = "_"))
  plot(spectrum)
  dev.off() 
}

#Raw data
compound = arrange(pivot_wider(rbindlist(exptlist, idcol = T), names_from = .id, values_from = Intensity),Mass)
colnames(compound) = CElist
compound$MassIonSum = rowSums(compound[,-1], na.rm=T)
compound = compound %>% arrange(desc(MassIonSum))
compound = rbind(compound, c(0,as.numeric(colSums(compound[,-1], na.rm=T))))

write.csv(compound,paste("MSMS", instrument[model+1], compoundID,"rawdata_fragmentation.csv",sep="_"),row.names = F)
compound_melted = compound[-nrow(compound),] %>% 
  mutate(mass = replace(mass, mass == "166.1"| mass == "166.15" | mass == "166.00", 166.05)) %>%
  group_by(mass) %>%
  summarise(mass, across(everything(), sum, na.rm=T)) %>%
  ungroup() %>%
  #filter(MassIonSum > max(MassIonSum,na.rm = T)/50) %>%
  filter(mass < 167) %>%
  arrange(desc(MassIonSum)) %>%
  select(-MassIonSum) %>%
  slice_head(n = 7) %>%
  melt(id="mass")
raw_frags = ggplot(data = compound_melted, aes(x = unfactor(variable), y = value, color = as.factor(mass))) +
  scale_color_npg() + scale_fill_npg() +
  geom_point() +  geom_line() +
  labs (x = "Collision energy (eV)", 
        y = "Intensity", 
        title = str_wrap(paste(instrument[model], "-", compoundID, "positive MSMS fragmentation (raw)", sep = " "), width = 60),
        color = "m/z")

png(filename=paste("MSMS", instrument[model+1], compoundID, "raw fragmentation.png", sep = " "))
plot(raw_frags)
dev.off() 

#Normalised
compound_normalised = cbind(compound$mass,compound[,-1]*100/c(compound[nrow(compound),-1]))
#if (model == 7) {compound_normalised = cbind(compound$mass,compound_normalised[,-1]*100/c(compound[nrow(compound),-1]))
#                 compound_normalised = cbind(compound_normalised[1],compound_normalised[-1]*54.0389)
#}
colnames(compound_normalised)[1] = "mass"

write.csv(compound_normalised,paste("MSMS", instrument[model+1], compoundID,"normalised_fragmentation.csv",sep="_"),row.names = F)
compound_normalised_melted = compound_normalised[-nrow(compound_normalised),] %>% 
  mutate(mass = replace(mass, mass == "166.1"| mass == "166.15" | mass == "166.00", 166.05)) %>%
  group_by(mass) %>%
  summarise(mass, across(everything(), sum, na.rm=T)) %>%
  ungroup() %>%
  #filter(MassIonSum > max(MassIonSum,na.rm = T)/50) %>%
  filter(mass < 167) %>%
  arrange(desc(MassIonSum)) %>%
  select(-MassIonSum) %>%
  slice_head(n=7) %>%
  melt(id="mass")

norm_frags = ggplot(data = compound_normalised_melted, aes(x = unfactor(variable), y = value, color = as.factor(mass))) +
  scale_color_npg() + scale_fill_npg() +
  geom_point() +  geom_line() +
  labs (x = "Collision energy (eV)", 
        y = "Intensity (%)", 
        title = str_wrap(paste(instrument[model], "-", compoundID, "positive MSMS fragmentation (normalised)", sep = " "), width = 60),
        color = "m/z")
png(filename=paste("MSMS", instrument[model+1], compoundID, "normalised fragmentation.png", sep = " "))
plot(norm_frags)
dev.off() 