@echo off
set nTurbines=9

for /L %%n in (1,1,%nTurbines%) do (
        start "" /B cmd /C D:\data\40_mahfouz\TurbSim\TurbSim\TurbSim.exe D:\data\40_mahfouz\FAST_Farm_template\Turbsim_Case\Cond00_v10.0_PL0.2_TI6\Case0_wdirp00\Seed_0\TurbSim\HighT%%n.inp > D:\data\40_mahfouz\FAST_Farm_template\Turbsim_Case\Cond00_v10.0_PL0.2_TI6\Case0_wdirp00\Seed_0\TurbSim\log.hight%%n.seed0.txt 2>&1
)

echo Script execution completed
