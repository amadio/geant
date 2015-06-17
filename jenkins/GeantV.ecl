-project=GeantV,getenv("PROJECT_ROOT")
-report_output={protobuf,getenv("PB_OUTPUT")}
-report_enabled={text,false},{protobuf,true}
-report_macros={protobuf,10}

-eval-file=GEANT_rules_config.ecl

-file_category={system,"^(/|<)"}

-locations={hide,{},system, "Disregard non-GEANT sources."}
-files={hide,system, "Disregard non-GEANT sources."}
-source_files={hide,system,"Disregard non-GEANT sources."}

-locations+={hide,{context},all, "Context locations are uninteresting."}
