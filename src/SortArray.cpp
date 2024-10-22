/*******************************************************************************
 * src/SortArray.cpp
 *
 * SortArray represents a simple array, which is displayed by WSortView to the
 * user is real-time.
 *
 *******************************************************************************
 * Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include "SortArray.h"
#include "SortAlgo.h"

#include <algorithm>
#include <random>
#include <cmath>

extern void SoundAccess(size_t i);

// *****************************************************************************
// *** Comparisons of ArrayItems

size_t g_compare_count = 0;

size_t g_access_count = 0;

size_t m_swaps = 0;

void ArrayItem::OnAccess(const ArrayItem& a)
{
    SoundAccess(a.get_direct());
}

void ArrayItem::OnComparison(const ArrayItem& a, const ArrayItem& b)
{
    ++g_compare_count;

    SoundAccess(a.get_direct());
    SoundAccess(b.get_direct());
}

// *****************************************************************************
// *** SortArray

SortArray::SortArray()
    : m_calc_inversions(false),
      m_delay(nullptr)
{
}

void SortArray::OnAlgoLaunch(const AlgoEntry& ae)
{
    if (size() <= ae.inversion_count_limit)
    {
        m_calc_inversions = true;
        RecalcInversions();
    }
    else
    {
        m_calc_inversions = false;
        m_inversions = -1;
    }
}

void SortArray::ResetArray(size_t size)
{
    m_array.resize(size, ArrayItem(0));
    m_mark.resize(size);
}

void SortArray::FinishFill()
{
    ASSERT(m_array.size() > 0);

    // calculate max value in array
    m_array_max = m_array[0].get_direct();
    for (size_t i = 1; i < size(); ++i)
    {
        if (m_array_max < m_array[i].get_direct())
            m_array_max = m_array[i].get_direct();
    }

    // reset access markers
    unmark_all();
    unwatch_all();

    // reset counters and info
    m_is_sorted = false;
    g_access_count = 0;
    g_compare_count = 0;
    m_calc_inversions = true;
    m_swaps = 0;

    RecalcInversions();
}

void SortArray::FillInputlist(wxArrayString& list)
{
    list.Add(_("Random Shuffle"));
    list.Add(_("Ascending"));
    list.Add(_("Descending"));
    list.Add(_("Near Sorted"));
    list.Add(_("25% Shuffled (Tail)"));
    list.Add(_("25% Shuffled (Head)"));
    list.Add(_("50% Shuffled (Tail)"));
    list.Add(_("50% Shuffled (Head)"));
    list.Add(_("75% Shuffled (Tail)"));
    list.Add(_("75% Shuffled (Head)"));
    list.Add(_("Shuffled Cubic"));
    list.Add(_("Shuffled Quintic"));
    list.Add(_("Shuffled n-2 Equal"));
    list.Add(_("Pipe Organ"));
    list.Add(_("Mirrored Organ"));
    list.Add(_("Wave"));
    list.Add(_("Sawtooth"));
    list.Add(_("Reverse Sawtooth"));
    list.Add(_("Many Similar"));
    list.Add(_("Quicksort Killer"));
    list.Add(_("Spike"));
    list.Add(_("Ribbon"));
    list.Add(_("Max Heapified"));
    list.Add(_("Flipped Min Heapified"));
}

void minheapify(std::vector<ArrayItem>& m_array, int n, int i) {
    int smallest = i, left = 2 * i + 1, right = 2 * i + 2;
    if (left < n && m_array[left] < m_array[smallest]) { smallest = left; }
    if (right < n && m_array[right] < m_array[smallest]) { smallest = right; }
    if (smallest != i) 
    {
        ArrayItem temp = m_array[i];
        m_array[i] = m_array[smallest];
        m_array[smallest] = temp;
        minheapify(m_array, n, smallest);
    }
}

void heapify(std::vector<ArrayItem>& m_array, int n, int i)
{
    int largest = i, left = 2 * i + 1, right = 2 * i + 2;
    if (left < n && m_array[left] > m_array[largest]) { largest = left; }
    if (right < n && m_array[right] > m_array[largest]) { largest = right; }
    if (largest != i)
    {
        ArrayItem temp = m_array[i];
        m_array[i] = m_array[largest];
        m_array[largest] = temp;
        heapify(m_array, n, largest);
    }
}

void SortArray::FillData(unsigned int schema, size_t arraysize)
{
    if (arraysize == 0) arraysize = 1;
    { ResetArray(arraysize); }
    std::random_device rd;
    std::mt19937 g(rd());
    switch (schema)
    {
        case 0:  // Shuffle of [1,n]
        {
            for (size_t i = 0; i < m_array.size(); ++i) { m_array[i] = ArrayItem(i + 1); }
            std::shuffle(m_array.begin(), m_array.end(), g);
            break;
        }
        case 1:  // Ascending [1,n]
        {
            for (size_t i = 0; i < m_array.size(); ++i) { m_array[i] = ArrayItem(i + 1); }
            break;
        }
        case 2:  // Descending
        {
            for (size_t i = 0; i < m_array.size(); ++i) { m_array[i] = ArrayItem(m_array.size() - i); }
            break;
        }
        case 3:  // Near Sorted
        {
            for (size_t i = 0; i < m_array.size(); ++i) { m_array[i] = ArrayItem(i + 1); }
            ArrayItem temp1 = m_array[0];
            m_array[0] = m_array[m_array.size() - 1];
            m_array[m_array.size() - 1] = temp1;
            break;
        }
        case 4:  // 25% Shuffled
        {
            std::vector<ArrayItem>::iterator it1 = m_array.begin();
            size_t half = m_array.size() / 2;
            for (size_t i = 0; i < m_array.size(); ++i)
            {
                m_array[i] = ArrayItem(i + 1);
                if (i <= half + (half / 2))
                { ++it1; }
            }
            std::shuffle(it1, m_array.end(), g);
            break;
        }
        case 5:  // 25% Sorted, Head
        {
            std::vector<ArrayItem>::iterator it = m_array.begin();
            for (size_t i = 0; i < m_array.size(); ++i)
            {
                m_array[i] = ArrayItem(i + 1);
                if (i <= (m_array.size() / 4))
                { ++it; }
            }
            std::shuffle(m_array.begin(), it, g);
            break;
        }
        case 6:  // 50% Sorted
        {
            std::vector<ArrayItem>::iterator it2 = m_array.begin();
            for (size_t i = 0; i < m_array.size(); ++i)
            {
                m_array[i] = ArrayItem(i + 1);
                if (i <= (m_array.size() / 2) - 1)
                { ++it2; }
            }
            std::shuffle(it2, m_array.end(), g);
            break;
        }
        case 7:  // 50% Sorted, Head
        {
            std::vector<ArrayItem>::iterator it3 = m_array.begin();
            for (size_t i = 0; i < m_array.size(); ++i)
            {
                m_array[i] = ArrayItem(i + 1);
                if (i <= (m_array.size() / 2) - 1)
                { ++it3; }
            }
            std::shuffle(m_array.begin(), it3, g);
            break;
        }
        case 8:  // 75% Shuffled
        {
            std::vector<ArrayItem>::iterator it4 = m_array.begin();
            for (size_t i = 0; i < m_array.size(); ++i)
            {
                m_array[i] = ArrayItem(i + 1);
                if (i <= (m_array.size() / 4))
                { ++it4; }
            }
            std::shuffle(it4, m_array.end(), g);
            break;
        }
        case 9:  // 75% Shuffled, Head
        {
            std::vector<ArrayItem>::iterator it5 = m_array.begin();
            size_t half = m_array.size() / 2;
            for (size_t i = 0; i < m_array.size(); ++i)
            {
                m_array[i] = ArrayItem(i + 1);
                if (i <= half + (half / 2))
                { ++it5; }
            }
            std::shuffle(m_array.begin(), it5, g);
            break;
        }
        case 10:  // Cubic skew of [1,n]
        {
            for (size_t i = 0; i < m_array.size(); ++i)
            {
                // normalize to [-1,+1]
                double x = (2.0 * (double)i / m_array.size()) - 1.0;
                // calculate x^3
                double v = x * x * x;
                // normalize to array size
                double w = (v + 1.0) / 2.0 * arraysize + 1;
                // decrease resolution for more equal values
                w /= 3.0;
                m_array[i] = ArrayItem(w + 1);
            }
            std::shuffle(m_array.begin(), m_array.end(), g);
            break;
        }
        case 11:  // Quintic skew of [1,n]
        {
            for (size_t i = 0; i < m_array.size(); ++i)
            {
                // normalize to [-1,+1]
                double x = (2.0 * (double)i / m_array.size()) - 1.0;
                // calculate x^5
                double v = x * x * x * x * x;
                // normalize to array size
                double w = (v + 1.0) / 2.0 * arraysize + 1;
                // decrease resolution for more equal values
                w /= 3.0;
                m_array[i] = ArrayItem(w + 1);
            }
            std::shuffle(m_array.begin(), m_array.end(), g);
            break;
        }
        case 12:  // shuffled n-2 equal values in [1,n]
        {
            m_array[0] = ArrayItem(1);
            for (size_t i = 1; i < m_array.size() - 1; ++i) { m_array[i] = ArrayItem(arraysize / 2 + 1); }
            m_array[m_array.size() - 1] = ArrayItem(arraysize);
            std::shuffle(m_array.begin(), m_array.end(), g);
            break;
        }
        case 13:  // Pipe organ (1, 1, 2, 2, 1, 1)
        {
            size_t n = m_array.size();
            if (n % 2 == 0)
            {
                int val = 1;
                for (size_t i = 0; i < n / 2; ++i)
                {
                    m_array[i] = ArrayItem(val); ++val;
                }
                val = n / 2;
                for (size_t i = n / 2; i <= n - 1; ++i)
                {
                    m_array[i] = ArrayItem(val); --val;
                }
            }
            else
            {
                int val = 1;
                for (size_t i = 0; i <= n / 2; ++i)
                {
                    m_array[i] = ArrayItem(val); ++val;
                }
                val = n / 2;
                for (size_t i = (n / 2) + 1; i <= n - 1; ++i)
                {
                    m_array[i] = ArrayItem(val); --val;
                }
            }
            break;
        }
        case 14:  // Mirrored organ (3, 2, 1, 1, 2, 3)
        {
            size_t n = m_array.size();
            if (n % 2 == 0)
            {
                int val = n / 2;
                for (size_t i = 0; i < n / 2; ++i)
                {
                    m_array[i] = ArrayItem(val); --val;
                }
                val = 1;
                for (size_t i = n / 2; i <= n - 1; ++i)
                {
                    m_array[i] = ArrayItem(val); ++val;
                }
            }
            else
            {
                int val = n / 2;
                for (size_t i = 0; i <= n / 2; ++i)
                {
                    m_array[i] = ArrayItem(val);
                    if (val >= 2) { --val; }
                }
                val = 1;
                for (size_t i = (n / 2) + 1; i <= n - 1; ++i)
                {
                    m_array[i] = ArrayItem(val); ++val;
                }
            }
            break;
        }
        case 15:  // Wave
        {
            double n = double(m_array.size());
            double pi = 3.14159265358979323846;
            for (size_t i = 0; i < m_array.size(); ++i)
            {
                double x = i / n * 3 * pi * 5;
                double sineVal = sin(x);
                int val = std::round((sineVal + 1) * 100);
                m_array[i] = ArrayItem(val + 1);
            }
            break;
        }
        case 16:  // Sawtooth
        {
            size_t n = m_array.size(), teeth;
            if (n % 5 == 0) { teeth = 5; }
            else if (n % 4 == 0) { teeth = 4; }
            else if (n % 3 == 0) { teeth = 3; }
            else { teeth = 2; }
            int max = n / teeth;
            int count = 1;
            for (size_t i = 0; i < n; ++i)
            {
                if (count > max) { count = 1; }
                m_array[i] = ArrayItem(count);
                ++count;
            }
            if (teeth == 2 && m_array[n - 1] == 1)
            {
                size_t m = n - 1;
                while (m_array[m - 1] > m_array[m])
                {
                    ArrayItem temp5 = m_array[m - 1];
                    m_array[m - 1] = m_array[m];
                    m_array[m] = temp5;
                    --m;
                }
            }
            break;
        }
        case 17:  // Reverse Sawtooth
        {
            size_t n = m_array.size(), teeth;
            if (n % 5 == 0) { teeth = 5; }
            else if (n % 4 == 0) { teeth = 4; }
            else if (n % 3 == 0) { teeth = 3; }
            else { teeth = 2; }
            int max = n / teeth;
            int count = max;
            for (size_t i = 0; i < n; ++i)
            {
                if (count <= 0) { count = max; }
                m_array[i] = ArrayItem(count);
                --count;
            }
            if (teeth == 2 && m_array[n - 1] == max)
            {
                size_t m = n - 1;
                while (m_array[m - 1] < m_array[m])
                {
                    ArrayItem temp4 = m_array[m - 1];
                    m_array[m - 1] = m_array[m];
                    m_array[m] = temp4;
                    --m;
                }
            }
            break;
        }
        case 18:  // Many Similar
        {
            size_t group_count = 0, size = m_array.size();
            if (size % 10 == 0) { group_count = 10; }
            else if (size % 9 == 0) { group_count = 9; }
            else if (size % 8 == 0) { group_count = 8; }
            else if (size % 7 == 0) { group_count = 7; }
            else if (size % 6 == 0) { group_count = 6; }
            else if (size % 5 == 0) { group_count = 5; }
            else if (size % 4 == 0) { group_count = 4; }
            else if (size % 2 != 0) { group_count = 3; }
            else { group_count = 2; }
            size_t n = m_array.size();
            if (n <= 10)
            {
                if (n % 2 != 0) { group_count = 3; }
                else { group_count = 2; }
            }
            size_t repeat = 1;
            int val = 1;
            for (size_t i = 0; i < size; ++i)
            {
                if (repeat > group_count)
                {
                    ++val; repeat = 1;
                }
                m_array[i] = ArrayItem(val);
                ++repeat;
            }
            std::shuffle(m_array.begin(), m_array.end(), g);
            break;
        }
        case 19:  // Quicksort Killer
        {
            int currentLen = m_array.size();
            for (int i = 0; i < currentLen; ++i)
            {
                m_array[i] = ArrayItem(i + 1);
            }
            for (int j = currentLen - currentLen % 2 - 2, i = j - 1; i >= 0; i -= 2, j--)
            {
                ArrayItem temp3 = m_array[i];
                m_array[i] = m_array[j];
                m_array[j] = temp3;
            }
            break;
        }
        case 20:  // Spike
        {
            size_t n = m_array.size();
            int spike, val = 1;
            if (n % 10 == 0) { spike = 5; }
            else if (n % 8 == 0) { spike = 4; }
            else if (n % 6 == 0) { spike = 3; }
            else { spike = 2; }
            int max = n / (spike * 2);
            if (max == 1) { max = n / spike; }
            for (size_t i = 0; i < n; ++i)
            {
                while (val <= max && i < n)
                {
                    m_array[i] = ArrayItem(val); ++val; ++i;
                }
                if (n % 2 == 0) { val = max; }
                else { val = max - 1; }
                while (val > 0 && i < n)
                {
                    m_array[i] = ArrayItem(val); --val;
                    if (val != 0) { ++i; }
                }
                val = 1;
            }
            size_t start, end;
            start = end = m_array.size() - 1;
            while (m_array[start - 1] < m_array[start])
            {
                --start;
            }
            for (; start <= end; ++start)
            {
                ArrayItem key = m_array[start];
                size_t j = start;
                while (j >= 1 && m_array[j - 1] <= key)
                {
                    m_array[j] = m_array[j - 1];
                    --j;
                }
                m_array[j] = key;
            }
            break;
        }
        case 21:  // Ribbon
        {
            int min = 1;
            int max = m_array.size();
            for (size_t i = 0; i < m_array.size(); ++i)
            {
                if (i % 2 == 0) { m_array[i] = ArrayItem(min); }
                else { m_array[i] = ArrayItem(max); }
                ++min; --max;
            }
            break;
        }
        case 22:  // Max Heapified
        {
            int n = m_array.size();
            for (int i = 0; i < n; ++i) { m_array[i] = ArrayItem(i + 1); }
            std::shuffle(m_array.begin(), m_array.end(), g);
            for (int i = n / 2 - 1; i >= 0; --i) { heapify(m_array, n, i); }
            break;
        }
        case 23:  // Flipped Min Heapified
        {
            int n = m_array.size();
            for (int i = 0; i < n; ++i) { m_array[i] = ArrayItem(i + 1); }
            std::shuffle(m_array.begin(), m_array.end(), g);
            for (int i = n / 2 - 1; i >= 0; --i) { minheapify(m_array, n, i); }
            std::reverse(m_array.begin(), m_array.end());
            break;
        }
        default:
            return FillData(0, arraysize);
            break;
    }
    FinishFill();
}

void SortArray::OnAccess()
{
    ++g_access_count;

    if (m_delay)
        m_delay->OnAccess();
}

bool SortArray::CheckSorted()
{
    unmark_all();
    // needed because iterator instrumentated algorithms may have changed the array
    RecalcInversions();
    mark(0);
    bool is_sorted = true;

    for (size_t i = 1; i < size(); ++i)
    {
        g_compare_count--; // dont count the following comparison
        if (get_nocount(i - 1) > get_nocount(i)) {
            wxLogError(_T("Result of sorting algorithm is incorrect!"));
            is_sorted = false;
            break;
        }
        mark(i);
    }

    unmark_all();
    return (m_is_sorted = is_sorted);
}

void SortArray::SetCalcInversions(bool on)
{
    // toggle boolean
    m_calc_inversions = on;

    if (!m_calc_inversions)
        m_inversions = -1;
}

void SortArray::ToggleCalcInversions()
{
    // toggle boolean
    SetCalcInversions(!m_calc_inversions);
}

void SortArray::RecalcInversions()
{
    if (!m_calc_inversions) {
        m_inversions = -1;
        return;
    }

    unsigned int inversions = 0;

    for (size_t i = 0; i < size(); ++i)
    {
        const ArrayItem& a = direct(i);

        for (size_t j = i+1; j < size(); ++j)
        {
            const ArrayItem& b = direct(j);

            if ( a.greater_direct(b) )
            {
                inversions++;
            }
        }
    }

    m_inversions = inversions;
}

void SortArray::AddInversions(size_t pos)
{
    if (!m_calc_inversions) {
        m_inversions = -1;
        return;
    }

    const ArrayItem& a = direct(pos);
    unsigned int inverses = 0;
    for (size_t i = 0; i < pos; ++i)
    {
        const ArrayItem& b = direct(i);
        if (b.greater_direct(a)) { ++inverses; }
    }

    m_inversions += inverses;
}

void SortArray::RemoveInversions(size_t pos)
{
    if (!m_calc_inversions) {
        m_inversions = -1;
        return;
    }

    const ArrayItem& a = direct(pos);
    unsigned int inverses = 0;
    for (size_t i = 0; i < pos; ++i)
    {
        const ArrayItem& b = direct(i);
        if (b.greater_direct(a)) { ++inverses; }
    }

    m_inversions -= inverses;
}

void SortArray::UpdateInversions(size_t i, size_t j)
{
    if (!m_calc_inversions) {
        m_inversions = -1;
        return;
    }
    if (m_inversions < 0) return RecalcInversions();

    if (i == j) return;

    unsigned int lo = i, hi = j;
    if (lo > hi) std::swap(lo, hi);

    const ArrayItem& ilo = m_array[lo];
    const ArrayItem& ihi = m_array[hi];
    int invdelta = 0;

    for (size_t k = lo + 1; k < hi; ++k)
    {
        if (m_array[k].less_direct(ilo))
            invdelta--;
        if (m_array[k].greater_direct(ilo))
            invdelta++;

        if (m_array[k].less_direct(ihi))
            invdelta++;
        if (m_array[k].greater_direct(ihi))
            invdelta--;
    }

    if (ilo.less_direct(ihi))
        invdelta++;
    if (ilo.greater_direct(ihi))
        invdelta--;

    m_inversions += invdelta;
}

size_t SortArray::GetRuns() const
{
    unsigned int runs = 1;

    for (size_t i = 1; i < size(); ++i)
    {
        const ArrayItem& a = direct(i-1);
        const ArrayItem& b = direct(i);

        if ( a.greater_direct(b) )
        {
            runs++;
        }
    }

    return runs;
}

short SortArray::InAccessList(ssize_t idx)
{
    if (idx < 0) return -1;

    signed color = -1;
    signed priority = -1;

    for (std::vector<Access>::iterator it = m_access_list.begin();
         it != m_access_list.end(); )
    {
        if (it->index != (size_t)idx) {
            ++it;
            continue;
        }

        if (it->priority >= priority)
        {
            priority = it->priority;
            color = it->color;
        }

        if (it->sustain == 0) {
            if (it->index == m_access1.index ||
                it->index == m_access2.index)
            {
                ++it;
            }
            else
            {
                it = m_access_list.erase(it);
            }
        }
        else {
            it->sustain--;
            ++it;
        }
    }

    return color;
}

unsigned short SortArray::InWatchList(ssize_t idx) const
{
    for (size_t i = 0; i < m_watch.size(); ++i)
    {
        if (m_watch[i].first == nullptr) continue;

        // compare watched value
        if (*m_watch[i].first != idx) continue;

        return m_watch[i].second;
    }
    return 0;
}

int SortArray::GetIndexColor(size_t idx)
{
    int clr, acl = InAccessList(idx);

    // select color
    if (idx == m_access1.index)
    {
        clr = m_access1.color;
    }
    else if (idx == m_access2.index)
    {
        clr = m_access2.color;
    }
    else if ( (clr = InWatchList(idx)) != 0 )
    {
        // clr already set
    }
    else if (m_mark[idx] != 0)
    {
        clr = m_mark[idx];
    }
    else if ( acl >= 0 )
    {
        clr = acl;
    }
    else
    {
        clr = 0;
    }

    return clr;
}

// *****************************************************************************
