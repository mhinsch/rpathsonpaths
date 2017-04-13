#ifndef PROPORTIONALPICK_H
#define PROPORTIONALPICK_H

#include <vector>

template<class FIT=double>
class ProportionalPick
	{
public:
	static const FIT & identity(const FIT & arg)
		{
		return arg;
		}

	ProportionalPick(const FIT & delta)
		: _delta(delta)
		{}
		
	template<class CONT>
	ProportionalPick(const FIT & delta, const CONT & cont)
		: _delta(delta)
		{
		setup(cont.begin(), cont.end());
		}

	template<class ITER, class FUNC>
	void setup(ITER start, ITER stop, FUNC fn)
		{
		FIT sum = FIT(0);
		_fitness.reserve(stop - start);
		
		for (ITER i=start; i!=stop; i++)
			{
			sum += fn(*i);
			_fitness.push_back(sum);
			}
		}

	template<class ITER>
	void setup(ITER start, ITER stop)
		{
		setup(start, stop, identity);
		}
	
	template<class RNG>
	size_t pick(RNG & rng) const
		{
		if (_fitness.back() <= _delta)
			return rng(_fitness.size());
		else
			{
			const FIT p = rng.outOf(FIT(0), _fitness.back());
			return lower_bound
				(_fitness.begin(), _fitness.end(), p) - _fitness.begin();
			}
		}
		
protected:
	std::vector<FIT> _fitness;
	// have to make it an instance variable
	// since floats can't be template parameters
	const FIT _delta;
	};

template<class FIT, class CONT, class RNG>
size_t propPick(const FIT & delta, const CONT & cont, RNG & rng)
	{
	ProportionalPick<FIT> pick(delta, cont);

	return pick.pick(rng);
	}

#endif	//PROPORTIONALPICK_H
