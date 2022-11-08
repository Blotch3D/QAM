using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework.Content;
using Blotch;
using System.Threading;
using System.Numerics;

namespace Blotch.Qam
{
    /// <summary>
    /// A quadrature amplitude modulator/demodulator. You can specify any constellation.
    /// Instantiate one of these on the transmitter to modulate symbols into a QAM signal suitable for sending over a channel.
    /// Instantiate one of these on the receiver (after acquisition and tracking) using the same parameters to demodulate
    /// the received signal into its symbols. 
    /// 
    /// The modulated carrier uses a signal sampling rate of four samples per full carrier wave. If your channel
    /// requires extra out-of-band rejection and you have the capability, you might want to upsample the signal before
    /// sending it.
    /// 
    /// Implement integration by upsampling and downsampling the states array in your app code.
    /// 
    /// A transmitter would typically add FEC to the source data, convert that channel code to symbols, and then convert
    /// the symbols to QAM signal by calling SymbolsToStates and then StatesToSignal.
    /// 
    /// A receiver must synchronize with symbol start, track the signal if the stream is long, perform AGC,
    /// demodulate by calling SignalToStates and then StatesToSymbols, and then decode the FEC. The process is designed so
    /// that you have access to the intermediate modulation states for doing these things optimally, integrating, etc.
    /// 
    /// </summary>
    public class QAM
    {
        public Dictionary<int, Complex> Constellation = null;
        /// <summary>
        /// See class description for details. This can be changed dynamically.
        /// </summary>
        int NumPoints;
        int NumSymbBits;
        List<int> Symbols;
        double SeverityRadius = 1;

        // Following are some typical constellations you can pass to the constructor. (You can also specify custom
        // constellations, or create one programmatically with CreateConstellation.

        public static Dictionary<int, Complex> QAM8_circular_7_1 = new Dictionary<int, Complex>()
            {
                { 3, new Complex(1, 0) },
                { 5, new Complex(0.6234898018587336, 0.7818314824680298) },
                { 4, new Complex(-0.22252093395631434, 0.9749279121818236) },
                { 6, new Complex(-0.900968867902419, 0.43388373911755823) },
                { 2, new Complex(-0.9009688679024191, -0.433883739117558) },
                { 0, new Complex(-0.2225209339563146, -0.9749279121818236) },
                { 1, new Complex(0.6234898018587334, -0.7818314824680299) },
                { 7, new Complex(0, 0) },
        };

        public static Dictionary<int, Complex> QAM16_rectangular = new Dictionary<int, Complex>()
            {
                { 15, new Complex(-3, 3) },
                { 7, new Complex(-1, 3) },
                { 5, new Complex(1, 3) },
                { 13, new Complex(3, 3) },
                { 14, new Complex(-3, 1) },
                { 6, new Complex(-1, 1) },
                { 4, new Complex(1, 1) },
                { 12, new Complex(3, 1) },
                { 10, new Complex(-3, -1) },
                { 2, new Complex(-1, -1) },
                { 0, new Complex(1, -1) },
                { 8, new Complex(3, -1) },
                { 11, new Complex(-3, -3) },
                { 3, new Complex(-1, -3) },
                { 1, new Complex(1, -3) },
                { 9, new Complex(3, -3) },
        };

        public static Dictionary<int, Complex> QAM32_rectangular = new Dictionary<int, Complex>()
            {
                { 21, new Complex(-3, 5) },
                { 17, new Complex(-1, 5) },
                { 1, new Complex(1, 5) },
                { 5, new Complex(3, 5) },
                { 7, new Complex(-5, 3) },
                { 23, new Complex(-3, 3) },
                { 19, new Complex(-1, 3) },
                { 3, new Complex(1, 3) },
                { 4, new Complex(3, 3) },
                { 20, new Complex(5, 3) },
                { 6, new Complex(-5, 1) },
                { 22, new Complex(-3, 1) },
                { 18, new Complex(-1, 1) },
                { 2, new Complex(1, 1) },
                { 0, new Complex(3, 1) },
                { 16, new Complex(5, 1) },
                { 14, new Complex(-5, -1) },
                { 30, new Complex(-3, -1) },
                { 26, new Complex(-1, -1) },
                { 10, new Complex(1, -1) },
                { 8, new Complex(3, -1) },
                { 24, new Complex(5, -1) },
                { 12, new Complex(-5, -3) },
                { 28, new Complex(-3, -3) },
                { 27, new Complex(-1, -3) },
                { 11, new Complex(1, -3) },
                { 9, new Complex(3, -3) },
                { 25, new Complex(5, -3) },
                { 29, new Complex(-3, -5) },
                { 31, new Complex(-1, -5) },
                { 15, new Complex(1, -5) },
                { 13, new Complex(3, -5) },
        };


        /// <summary>
        /// Creates a QAM modulator/demodulator. See class description for details.
        /// 
        /// </summary>
        /// <param name="constellation">A Dictionary that relates each symbol (Key) to its corresponding constellation
        /// point (Value). A constellation point is a (Real amplitude,Imaginary amplitude) pair.</param>
        public QAM(Dictionary<int, Complex> constellation)
        {
            if (constellation == null)
            {
                constellation = QAM32_rectangular;
            }
            Constellation = constellation;
            SetConstellation();
        }

        /// <summary>
        /// Shows the constellation in a 3D window
        /// </summary>

        Win3d Win = null;
        public void Visualize()
        {
            if (Win == null)
            {
                Win = Win3d.CreateWindow();
            }

            // Remove any old subsprites
            Win.EnqueueCommand((w) =>
            {
                while (Win.TopSprite.Count > 0)
                {
                    var s = Win.TopSprite.First().Value;
                    s.Dispose();
                    Win.TopSprite.Remove(s.Name);
                }
            });

            // add constellation points
            Win.EnqueueCommand((w) =>
            {
                var num = 0;
                foreach (var pnt in Constellation)
                {
                    var s = new BlSprite(Win.Graphics, $"p{num++}");
                    s.Text = $"0x{pnt.Key:X}";
                    s.TextFont = Win.Font;

                    s.LODs.Add(Win.SphereModel);

                    Win.TopSprite.Add(s);

                    s.Matrix = Matrix.CreateScale((float).04, (float).04, (float).04);
                    s.Matrix *= Matrix.CreateTranslation((float)pnt.Value.Real, (float)pnt.Value.Imaginary, 0);
                }
            });
        }

        /// <summary>
        /// Creates the constellation from the 'Constellation' class member
        /// </summary>
        void SetConstellation()
        {
            NumPoints = Constellation.Count;
            NumSymbBits = (int)(Math.Log10(NumPoints) / Math.Log10(2) + 1);
            Symbols = Constellation.Keys.ToList();
            var avgPower = CalcAveragePower();
            Scale(1/avgPower);
        }

        /// <summary>
        /// Call this if you want to scale the magnitudes of the constellation points
        /// </summary>
        /// <param name="factor"></param>
        public void Scale(double factor)
        {
            var c = Constellation;
            c = new Dictionary<int, Complex>();
            for (int n = 0; n < Constellation.Count; n++)
            {
                var p = Constellation[n];
                c.Add(n, p * factor);
            }
            Constellation = c;
        }

        /// <summary>
        /// Figures the average xmit/rcv power assuming each point has the same hit rate
        /// </summary>
        /// <returns></returns>
        public double CalcAveragePower()
        {
            var power = 0.0;
            foreach(var p in Constellation)
            {
                power += p.Value.Real * p.Value.Real + p.Value.Imaginary * p.Value.Imaginary;
            }
            return power / Constellation.Count;
        }

        public List<Complex> UpsampleStates(List<Complex> states, int count = 2)
        {
            var upStates = new List<Complex>();

            foreach (var s in states)
            {
                for(int i = 0; i < count; i++)
                {
                    upStates.Add(s);
                }
            }

            return upStates;
        }
        public List<Complex> DownsampleStates(List<Complex> states, int count = 2)
        {
            var downStates = new List<Complex>();

            int n = 0;
            for (int i = 0; i < states.Count / count; i++)
            {
                var c = new Complex();

                for(int j = 0; j < count; j++)
                {
                    c += states[n++];
                }
                c /= count;

                downStates.Add(c);
            }

            foreach (var s in states)
            {
                downStates.Add(s);
                downStates.Add(s);
            }

            return downStates;
        }


        /// <summary>
        /// Creates a rectangular or triangular constellation given how many points. Symbol assignments are not
        /// optimized. Use OptimizeSymbolLayout for optimal symbols.
        /// </summary>
        /// <param name="numPoints">How many points should be in the constellation</param>
        public void CreateConstellation(int numPoints, bool isTriangular = true, double offsetX = 0, double offsetY = 0)
        {
            // We create the triangular grid by first creating a much larger grid, and then trimmin git down until it
            // has the correct number of points

            var dim = (int)Math.Sqrt(numPoints) + 2;

            var vertSeparation = 1.0;

            if (isTriangular)
            {
                vertSeparation = Math.Sin(Math.PI / 3);
            }

            int pntNum = 0;

            // create the oversized triangular grid
            var points = new List<Complex>();
            for (var y = -dim; y < dim; y++)
            {
                var horizOffset = 0.0;
                if (isTriangular)
                {
                    horizOffset = .5 * (y & 1);
                }

                for (var x = (int)-dim; x < dim; x++)
                {
                    points.Add(new Complex(x + horizOffset + offsetX, vertSeparation * y + offsetY));
                }
            }

            // Sort by distance from origin.
            // This is slightly less efficient because we re-calculate distance squared every time, but in the big
            // picture that's OK
            points.Sort((a, b) =>
            {
                var aDistSqr = a.Real * a.Real + a.Imaginary * a.Imaginary;
                var bDistSqr = b.Real * b.Real + b.Imaginary * b.Imaginary;
                if (aDistSqr < bDistSqr) return -1;
                if (aDistSqr > bDistSqr) return 1;
                return 0;
            });

            // Now we eliminate all points except the closest numPoints
            points = points.GetRange(0, numPoints);

            Constellation.Clear();
            int n = 0;
            foreach (var pnt in points)
            {
                Constellation.Add(n++, pnt);
            }

            SetConstellation();
        }

        /// <summary>
        /// improve efficiency by re-creating Constellation each call
        /// </summary>
        /// <param name="numSeeds"></param>
        /// <param name="numFinalIterations"></param>
        /// <param name="severityRadius"></param>
        public void OptimizeSymbolLayout(int numSeeds = 50, int numFinalIterations = 50, double severityRadius = 2)
        {
            var bestFlipsSeverity = 1e200;
            var newConstellation = new Dictionary<int, Complex>();

            for (int seed = 0; seed < numSeeds; seed++)
            {
                var flipsSeverity = OptimizeSymbolLayoutOnce(seed, numFinalIterations, severityRadius);

                var maxFlipsSeverity = 0.0;
                for (int n = 0; n < flipsSeverity.Count; n++)
                {
                    var thisSeverity = flipsSeverity[n];
                    if (maxFlipsSeverity < thisSeverity)
                    {
                        maxFlipsSeverity = thisSeverity;
                    }
                }

                if (bestFlipsSeverity > maxFlipsSeverity)
                {
                    bestFlipsSeverity = maxFlipsSeverity;
                    newConstellation.Clear();
                    foreach(var pnt in Constellation)
                    {
                        newConstellation.Add(pnt.Key, new Complex(pnt.Value.Real, pnt.Value.Imaginary));
                    }
                }
            }
            Constellation = newConstellation;
        }


        /// <summary>
        /// Tries to optimize which symbol is associated with which constellation point, such that the fewest bits
        /// differ between adjacent symbols in the constellation, and more so the closer they are. This does NOT alter
        /// actual point positions.
        /// 
        /// Since symbol layouts have many metastable 'optimized' states, calling this only once with a single random
        /// seed will almost certainly not create the most optimal layout. Instead, call it many times with different
        /// random seeds and use the best result. You decide whether the result is best according to the returned
        /// report.
        /// 
        /// Optimization is performed by repeatedly swapping symbols between points when the swap would decrease the
        /// severity level between adjacent symbols.
        /// </summary>
        /// <param name="rndSeed">Seeds the random number generator, which is used to initially scramble the symbols,
        /// and then to pick symbol pairs to test whether they should be swapped.</param>
        /// <param name="numFinalIterations">Optimization will be considered complete when the constellation is iterated this
        /// many times without any symbols being swapped.</param>
        /// <param name="severityRadius">Equates to the noise level so that the severity can be associated with it. This
        /// is the width of the Gaussian noise curve. the larger this value, the wider the 'neighborhood' of adjacent
        /// symbols that should be checked for bit differences.</param>
        /// <returns>A list where an element index is the point index in Constellation.Keys.ToList(), and the element
        /// value is a severity indicator of the symbol bits that differ with adjacent (nearby) points and depending on
        /// how far away they are. From this you can find the average and peak bit differences to better judge the
        /// result.</returns>
        public List<double> OptimizeSymbolLayoutOnce(
            int rndSeed = 0,
            int numFinalIterations = 50,
            double severityRadius = 1
            )
        {
            SeverityRadius = severityRadius;
            numFinalIterations *= (Constellation.Count * Constellation.Count * Constellation.Count)/100;

           // Various initializations
           var rnd = new Random(rndSeed);
            List<Complex> Points;
            Points = Constellation.Values.ToList();
            var flippedBitsSeverity = new List<double>();

             // For each point we create a list of its nearby points
            var nearPnts = GetNearbyPointInfoList();

            // randomize constellation
            for (int m = 0; m < NumPoints; m++)
            {
                var p1 = rnd.Next(0, NumPoints);
                var p2 = rnd.Next(0, NumPoints);

                var tmp = Symbols[p1];
                Symbols[p1] = Symbols[p2];
                Symbols[p2] = tmp;
            }

            var totalSeverity = 0.0;

            // load list of flippedBitsSeverity (each element represents a single constellation point, and its value is
            // how many a severity related to total flipped bits and separation to with that point's neigbors
            for (int m = 0; m < NumPoints; m++)
            {
                var severity = GetFlippedBitsSeverity(Symbols[m], nearPnts[m]);
                flippedBitsSeverity.Add(severity);
                totalSeverity += severity;
            }

            var nonswapIterCnt = 0;

            // Repeatedly find pairs of points that can be swapped in order to lower bit differences between their
            // neigbors
            for (int n = 0; n < numFinalIterations; n++)
            {
                if (nonswapIterCnt > numFinalIterations) break;

                // indices of symbols to possibly swap
                var p1 = rnd.Next(0, NumPoints);
                var p2 = rnd.Next(0, NumPoints);

                // if they aren't equal
                if (p1 != p2)
                {
                    //
                    // get flip severity for all scenarios
                    //

                    // current flipped bits severity for p1
                    var q1 = flippedBitsSeverity[p1];

                    // current flipped bits severity for p2
                    var q2 = flippedBitsSeverity[p2];

                    // temporarily swap them so we can test it
                    var tmp = Symbols[p1];
                    Symbols[p1] = Symbols[p2];
                    Symbols[p2] = tmp;

                    // potential flipped bits severity for p1 if we swap
                    var t1 = GetFlippedBitsSeverity(Symbols[p1], nearPnts[p1]);

                    // potential flipped bits severity for p2 if we swap
                    var t2 = GetFlippedBitsSeverity(Symbols[p2], nearPnts[p2]);

                    // Is the swap bad (i.e. we swap them back)?
                    if (q1 + q2 < t1 + t2)
                    {
                        tmp = Symbols[p1];
                        Symbols[p1] = Symbols[p2];
                        Symbols[p2] = tmp;
                        nonswapIterCnt++;
                    }
                    else
                    {
                        nonswapIterCnt = 0;
                        totalSeverity = 0;
                        for (int m = 0; m < NumPoints; m++)
                        {
                            var flipsSeverity = GetFlippedBitsSeverity(Symbols[m], nearPnts[m]);
                            flippedBitsSeverity[m] = flipsSeverity;
                            totalSeverity += flipsSeverity;
                        }
                    }
                }
            }

            // Create a new constellation with the new symbol layout
            var optConst = new Dictionary<int, Complex>();
            for (int n = 0; n < NumPoints; n++)
            {
                optConst.Add(Symbols[n], Points[n]);
            }

            Constellation = optConst;

            return flippedBitsSeverity;
        }


        double GaussianDenom = Math.Sqrt(2 * Math.PI);
        /// <summary>
        /// Figures the total flipped bits severity between a given symbol and the several other nearby symbols.
        /// This is used by the optimizer
        /// </summary>
        /// <param name="symb">The symbol value of the constellation point in question</param>
        /// <param name="nearPnts">The distances (Values) to several other nearby symbols (Keys)</param>
        /// <returns>A severity factor. Higher means more bit differences between neighbors, and its worse for closer
        /// neighbors</returns>
        double GetFlippedBitsSeverity(int symb, Dictionary<int, double> nearPnts)
        {
            var flippedBitsSeverity = 0.0;
            foreach(var neighbor in nearPnts)
            {
                var x = neighbor.Value / SeverityRadius;
                var severity = Math.Pow(Math.E, -(x * x)/2)/GaussianDenom;
                var bitsFlipped = symb ^ Symbols[neighbor.Key];
                for(int m =0; m < NumSymbBits; m++)
                {
                    if((bitsFlipped & 1) == 1)
                    {
                        flippedBitsSeverity += severity;
                    }
                    bitsFlipped >>= 1;
                    if (bitsFlipped == 0) break;
                }

            }
            return flippedBitsSeverity;
        }

        /// <summary>
        /// used by OptimizeSymbolLayout
        /// </summary>
        /// <param name="constellation"></param>
        /// <param name="radius"></param>
        /// <returns></returns>
        List<Dictionary<int, double>> GetNearbyPointInfoList(
            Dictionary<int, Complex> constellation = null, double radius = 2)
        {
            if(constellation == null)
            {
                constellation = Constellation;
            }

            var numPnts = constellation.Count;
            var points = constellation.Values.ToList();

            // The key is a point, the value is a dictionary of the nearby points and their effect assuming Gaussian
            // noise
            var distances = new List<Dictionary<int, double>>();
            for (int n = 0; n < numPnts; n++)
            {
                // find distance to all other points
                var dists = new Dictionary<int, double>();
                for (int m = 0; m < numPnts; m++)
                {
                    // if its not the same point
                    if (m != n)
                    {
                        // calcl separation distance
                        var dx = points[n].Real - points[m].Real;
                        var dy = points[n].Imaginary - points[m].Imaginary;
                        var difSqr = dx * dx + dy * dy;

                        // add it to the list
                        dists.Add(m, difSqr);
                    }
                }

                // find nearest distance
                var nearestDist = 1e100;
                foreach(var d in dists)
                {
                    if (d.Value < nearestDist)
                    {
                        nearestDist = d.Value;
                    }
                }

                // now collect only the near points
                var nearPnts = new Dictionary<int, double>();
                foreach (var d in dists)
                {
                    if (nearestDist * radius >= d.Value)
                    {
                        nearPnts.Add(d.Key, d.Value);
                    }
                }

                // add it to the main dictionary
                distances.Add(nearPnts);
            }
            return distances;
        }


        /// <summary>
        /// Given a list of QAM states (four states per carrier wave), convert them to a list of symbols using the
        /// current constellation. The first state must be the first sample of a state (i.e. this
        /// does not do any tracking).
        /// </summary>
        /// <param name="states">A list of consecutive states, where each state represents 1/4 of a carrier wave and is
        /// a (real amplitude, imaginary amplitude) pair.</param>
        /// <returns>A tuple, where the first item is the list of demodulated symbols, the second item is an average of
        /// the square of each state's deviation from the closest constellation point, and the third item is the square
        /// of the largest deviation from its closest constellation point.</returns>
        public (List<int>, double, double) StatesToSymbols(List<Complex> states, int start = 0, int end = -1)
        {
            if (end < 0) end = states.Count;

            var symbols = new List<int>();

            var avgDeviationSqr = 0.0;
            var maxDeviationSqr = 0.0;

            for(int n = start;n<end;n++)
            {
                var state = states[n];

                var bestDifSqr = 1e100;
                KeyValuePair<int, Complex> bestPoint = new KeyValuePair<int, Complex>();

                foreach(var point in Constellation)
                {
                    var dx = state.Real - point.Value.Real;
                    var dy = state.Imaginary - point.Value.Imaginary;

                    var difSqr = dx * dx + dy * dy;

                    if(bestDifSqr > difSqr)
                    {
                        bestDifSqr = difSqr;
                        bestPoint = point;
                    }
                }

                if (maxDeviationSqr < bestDifSqr)
                {
                    maxDeviationSqr = bestDifSqr;
                }

                avgDeviationSqr += bestDifSqr;

                symbols.Add(bestPoint.Key);
            }

            return (symbols, avgDeviationSqr / symbols.Count, maxDeviationSqr);
        }

        /// <summary>
        /// Convert a list of data to a list of symbols
        /// </summary>
        /// <param name="data"></param>
        /// <returns></returns>
        public List<int> DataToSymbols(List<int> data)
        {
            var symbols = new List<int>();
            foreach (var dat in data)
            {
                symbols.Add(Symbols[dat]);
            }
            return symbols;
        }

        /// <summary>
        /// Convert a list of symbols to a list of data
        /// This may be inefficient
        /// </summary>
        /// <param name="symbols"></param>
        /// <returns></returns>
        public List<int> SymbolsToData(List<int> symbols)
        {
            var data = new List<int>();
            foreach (var dat in symbols)
            {
                data.Add(Symbols.IndexOf(dat));
            }
            return data;
        }


        /// <summary>
        /// Produces a list of QAM states from a sequence of symbols.
        ///
        /// 'states' are a list of sequential real+imaginary amplitudes, one for each symbol in the symbol stream.
        /// To produce a signal suitable for transmission, pass the states to the StatesToSignal method.
        /// </summary>
        /// <param name="symbols">Data must first be encoded into symbols with the DataToSymbols method, then those
        /// symbols are passed as this parameter. No symbol can be outside the bounds of the constellation.</param>
        /// <returns>Returns a list of (real, imaginary) pairs, one for each input symbol</returns>
        public List<Complex> SymbolsToStates(List<int> symbols)
        {
            return symbols.Select(s => Constellation[s]).ToList();
        }

        /// <summary>
        /// Converts a sequence of consecutive atomic QAM states into a signal suitable for transmission.
        /// </summary>
        /// <param name="states">A sequence of QAM states, where each state defines the real and imaginary
        /// amplitudes of a complex carrier wave.</param>
        /// <returns>A signal suitable for transmission, where there are four samples per wave.</returns>
        public List<double> StatesToSignal(List<Complex> states, int startIdx = 0, int endIdx = -1, List<double> signalBuf = null)
        {
            if(endIdx < 0)
            {
                endIdx = states.Count - 1;
            }

            var numSignalSamples = (endIdx + 1 - startIdx) * 4;

            var signal = signalBuf;
            if(signal == null)
            {
                signal = new List<double>();
            }

            for (int stateIdx=startIdx; stateIdx<=endIdx; stateIdx++)
            {
                var state = states[stateIdx];

                signal.Add(state.Real);
                signal.Add(state.Imaginary);
                signal.Add(-state.Real);
                signal.Add(-state.Imaginary);
            }

            return signal;
        }

        /// <summary>
        /// Converts a signal to a list of QAM states.
        /// </summary>
        /// <param name="signal">A sequence of signal levels, where each level is 1/4 of a full wave</param>
        /// <param name="gain">A multiplier for the signal</param>
        /// <param name="start">Offset of first state</param>
        /// <param name="end">Offset of last state plus one</param>
        /// <returns>A list of consecutive states, where each state is the (real, imaginary) pair for the given quarter
        /// wav</returns>
        public List<Complex> SignalToStates(List<double> signal, double gain = 1, int startIdx = 0, int endIdx = -1)
        {
            var signalIdx = startIdx;

            if (endIdx < 0)
            {
                endIdx = signal.Count() - 1;
            }

            // convert signal indices to states indices
            startIdx /= 4;
            endIdx /= 4;

            var states = new List<Complex>();


            for (int stateIdx = startIdx; stateIdx <= endIdx; stateIdx++)
            {
                var state = new Complex((signal[signalIdx] - signal[signalIdx + 2]) / 2, (signal[signalIdx + 1] - signal[signalIdx + 3]) / 2);

                states.Add(state);

                signalIdx += 4;
            }

            return states;
        }
    }
}









