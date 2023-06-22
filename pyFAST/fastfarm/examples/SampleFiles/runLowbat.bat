@echo off
set nSeeds=1

for /L %%n in (0,1,%nSeeds%) do (
    if %%n LSS %nSeeds% (
    start ""/B cmd /C D:\data\40_mahfouz\TurbSim\TurbSim\TurbSim.exe D:\data\40_mahfouz\FAST_Farm_template\Turbsim_Case\Cond00_v10.0_PL0.2_TI6\Seed_%%n\Low.inp > D:\data\40_mahfouz\FAST_Farm_template\Turbsim_Case\Cond00_v10.0_PL0.2_TI6\Seed_%%n\log.low.seed%%n.txt 2>&1 
) 

)

echo Script execution completed
