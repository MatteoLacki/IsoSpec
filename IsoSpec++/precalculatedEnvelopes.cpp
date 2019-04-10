/*
 *   Copyright (C) 2015-2019 Mateusz Łącki and Michał Startek.
 *
 *   This file is part of IsoSpec.
 *
 *   IsoSpec is free software: you can redistribute it and/or modify
 *   it under the terms of the Simplified ("2-clause") BSD licence.
 *
 *   IsoSpec is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *   You should have received a copy of the Simplified BSD Licence
 *   along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
 */

#include "precalculatedEnvelopes.h"

namespace IsoSpec
{

template<typename T> void call_init(T* tabulator, Iso&& iso, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs)
{
    if(tgetlProbs)
    {
        if(tgetMasses)
        {
            if(tgetProbs)
            {
                if(tgetConfs)
                    tabulator->template init<true, true, true, true>(std::move(iso));
                else
                    tabulator->template init<true, true, true, false>(std::move(iso));
            }
            else
            {
                if(tgetConfs)
                    tabulator->template init<true, true, false, true>(std::move(iso));
                else
                    tabulator->template init<true, true, false, false>(std::move(iso));
            }
        }
        else
        {
            if(tgetProbs)
            {
                if(tgetConfs)
                    tabulator->template init<true, false, true, true>(std::move(iso));
                else
                    tabulator->template init<true, false, true, false>(std::move(iso));
            }
            else
            {
                if(tgetConfs)
                    tabulator->template init<true, false, false, true>(std::move(iso));
                else
                    tabulator->template init<true, false, false, false>(std::move(iso));
            }
        }
    }
    else
    {
        if(tgetMasses)
        {
            if(tgetProbs)
            {
                if(tgetConfs)
                    tabulator->template init<false, true, true, true>(std::move(iso));
                else
                    tabulator->template init<false, true, true, false>(std::move(iso));
            }
            else
            {
                if(tgetConfs)
                    tabulator->template init<false, true, false, true>(std::move(iso));
                else
                    tabulator->template init<false, true, false, false>(std::move(iso));
            }
        }
        else
        {
            if(tgetProbs)
            {
                if(tgetConfs)
                    tabulator->template init<false, false, true, true>(std::move(iso));
                else
                    tabulator->template init<false, false, true, false>(std::move(iso));
            }
            else
            {
                if(tgetConfs)
                    tabulator->template init<false, false, false, true>(std::move(iso));
                else
                    tabulator->template init<false, false, false, false>(std::move(iso));
            }
        }
    }
}

template void call_init<LayeredTabulator>(LayeredTabulator* tabulator, Iso&& iso, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs);
template void call_init<ThresholdTabulator>(ThresholdTabulator* tabulator, Iso&& iso, bool tgetlProbs, bool tgetMasses, bool tgetProbs, bool tgetConfs);


} // namespace IsoSpec
