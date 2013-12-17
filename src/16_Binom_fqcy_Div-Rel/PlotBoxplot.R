print(sample_size3 <- table(GLMtab_all2$Family, interaction(GLMtab_all2$trimmed, GLMtab_all2$CpDup, GLMtab_all2$Phylog_lvl)))

boxplot(GLMtab_all2$Fqcy_all~interaction(GLMtab_all2$Family, GLMtab_all2$trimmed, GLMtab_all2$CpDup, GLMtab_all2$Phylog_lvl))
levels(interaction(GLMtab_all2$Family, GLMtab_all2$trimmed, GLMtab_all2$CpDup, GLMtab_all2$Phylog_lvl))
#~ [1] "Control.No.FALSE.divergent"  "Gr.No.FALSE.divergent"      
#~ [3] "Or.No.FALSE.divergent"       "P450.No.FALSE.divergent"    
#~ [5] "Control.Yes.FALSE.divergent" "Gr.Yes.FALSE.divergent"     
#~ [7] "Or.Yes.FALSE.divergent"      "P450.Yes.FALSE.divergent"   
#~ [9] "Control.No.TRUE.divergent"   "Gr.No.TRUE.divergent"       
#~[11] "Or.No.TRUE.divergent"        "P450.No.TRUE.divergent"     
#~[13] "Control.Yes.TRUE.divergent"  "Gr.Yes.TRUE.divergent"      
#~[15] "Or.Yes.TRUE.divergent"       "P450.Yes.TRUE.divergent"    
#~[17] "Control.No.FALSE.related"    "Gr.No.FALSE.related"        
#~[19] "Or.No.FALSE.related"         "P450.No.FALSE.related"      
#~[21] "Control.Yes.FALSE.related"   "Gr.Yes.FALSE.related"       
#~[23] "Or.Yes.FALSE.related"        "P450.Yes.FALSE.related"     
#~[25] "Control.No.TRUE.related"     "Gr.No.TRUE.related"         
#~[27] "Or.No.TRUE.related"          "P450.No.TRUE.related"       
#~[29] "Control.Yes.TRUE.related"    "Gr.Yes.TRUE.related"        
#~[31] "Or.Yes.TRUE.related"         "P450.Yes.TRUE.related"      
boxplot(GLMtab_all2$Fqcy_all~interaction(GLMtab_all2$Family, GLMtab_all2$trimmed, GLMtab_all2$CpDup, GLMtab_all2$Phylog_lvl), col=colo)
abline(v=16.5, lty=1)
abline(v=24.5, lty=2)
abline(v=8.5, lty=2)
abline(v=4.5, lty=3)
abline(v=12.5, lty=3)
abline(v=20.5, lty=3)
abline(v=28.5, lty=3)
