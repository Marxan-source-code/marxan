#pragma once

#include <string>

namespace marxan {
    using namespace std;

    typedef struct srunoptions
    {
        int CalcPenaltiesOn;
        int HeuristicOn;
        int ThermalAnnealingOn;
        int HillClimbingOn;
        int ItImpOn;
        int TwoStepHillClimbingOn;
        
        string algorithm_description()
        {
            string res;
            if(ThermalAnnealingOn)
                res+= "Thermal Annealing";
            if(HeuristicOn)
            {
                if(!res.empty())
                    res+=", ";
                res+= "Heuristic";
            }
            if(HillClimbingOn)
            {
                if(!res.empty())
                    res+=", ";
                res+= "Hill Climbing";
            }
            if(TwoStepHillClimbingOn)
            {
                if(!res.empty())
                    res+=", ";
                res+= "Two Step Hill Climbing";
            }
            if(ItImpOn)
            {
                if(!res.empty())
                    res+=", ";
                res+= "Iterative Improvement";
            }
            //repace last comma with "and" 
            auto pos = res.find_last_of(',');
            if(pos != string::npos)
                res.replace(pos, 1, " and");
            return res;
        }

        void setDefaultRunOptions(int runopts)
        {
            switch (runopts)
            {
            case 0:
                CalcPenaltiesOn = 1;
                ThermalAnnealingOn = 1;
                HillClimbingOn = 0;
                HeuristicOn = 1;
                ItImpOn = 0;
                TwoStepHillClimbingOn = 0;
                break;
            case 1:
                CalcPenaltiesOn = 1;
                ThermalAnnealingOn = 1;
                HillClimbingOn = 0;
                HeuristicOn = 0;
                ItImpOn = 1;
                TwoStepHillClimbingOn = 0;
                break;
            case 2:
                CalcPenaltiesOn = 1;
                ThermalAnnealingOn = 1;
                HillClimbingOn = 0;
                HeuristicOn = 1;
                ItImpOn = 1;
                TwoStepHillClimbingOn = 0;
                break;
            case 3:
                CalcPenaltiesOn = 0;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 0;
                HeuristicOn = 1;
                ItImpOn = 0;
                TwoStepHillClimbingOn = 0;
                break;
            case 4:
                CalcPenaltiesOn = 0;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 0;
                HeuristicOn = 0;
                ItImpOn = 1;
                TwoStepHillClimbingOn = 0;
                break;
            case 5:
                CalcPenaltiesOn = 0;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 0;
                HeuristicOn = 1;
                ItImpOn = 1;
                TwoStepHillClimbingOn = 0;
                break;
            case 6:
                CalcPenaltiesOn = 1;
                ThermalAnnealingOn = 1;
                HillClimbingOn = 0;
                HeuristicOn = 0;
                ItImpOn = 0;
                TwoStepHillClimbingOn = 0;
                break;
            case 7:
                CalcPenaltiesOn = 1;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 0;
                HeuristicOn = 0;
                ItImpOn = 0;
                TwoStepHillClimbingOn = 0;
                break;
            case 8:
                CalcPenaltiesOn = 0;
                ThermalAnnealingOn = 1;
                HillClimbingOn = 0;
                HeuristicOn = 0;
                ItImpOn = 0;
                TwoStepHillClimbingOn = 0;
                break;
            case 9:
                CalcPenaltiesOn = 0;
                ThermalAnnealingOn = 1;
                HillClimbingOn = 0;
                HeuristicOn = 0;
                ItImpOn = 1;
                TwoStepHillClimbingOn = 0;
                break;
            case 10:
                CalcPenaltiesOn = 0;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 0;
                HeuristicOn = 0;
                ItImpOn = 1;
                TwoStepHillClimbingOn = 0;
                break;
            case 11:
                CalcPenaltiesOn = 1;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 1;
                HeuristicOn = 0;
                ItImpOn = 0;
                TwoStepHillClimbingOn = 0;
                break;
            case 12:
                CalcPenaltiesOn = 1;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 1;
                HeuristicOn = 0;
                ItImpOn = 1;
                TwoStepHillClimbingOn = 0;
                break;
            case 13:
                CalcPenaltiesOn = 0;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 1;
                HeuristicOn = 0;
                ItImpOn = 0;
                TwoStepHillClimbingOn = 0;
                break;
            case 14:
                CalcPenaltiesOn = 0;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 1;
                HeuristicOn = 0;
                ItImpOn = 1;
                TwoStepHillClimbingOn = 0;
                break;
            case 15:
                CalcPenaltiesOn = 1;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 0;
                HeuristicOn = 0;
                ItImpOn = 1;
                TwoStepHillClimbingOn = 0;
                break;
            case 16:
                CalcPenaltiesOn = 1;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 0;
                HeuristicOn = 0;
                ItImpOn = 0;
                TwoStepHillClimbingOn = 1;
                break;
            case 17:
                CalcPenaltiesOn = 1;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 1;
                HeuristicOn = 0;
                ItImpOn = 0;
                TwoStepHillClimbingOn = 1;
                break;
            case 18:
                CalcPenaltiesOn = 1;
                ThermalAnnealingOn = 1;
                HillClimbingOn = 1;
                HeuristicOn = 0;
                ItImpOn = 0;
                TwoStepHillClimbingOn = 1;
                break;

            default:
                CalcPenaltiesOn = 0;
                ThermalAnnealingOn = 0;
                HillClimbingOn = 0;
                HeuristicOn = 0;
                ItImpOn = 0;
                TwoStepHillClimbingOn = 0;
                break;
            }
        } // setDefaultRunOptions
        

    } srunoptions;

    // set default run options based on the selection algorithm chosen



} // namespace marxan