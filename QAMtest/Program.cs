using Blotch.Qam;
using System.Numerics;

namespace QAMtest
{
    class Program
    {
        /// used to generate random src data
        static Random Rnd = new Random(0);

        static void Main(string[] args)
        {
            TestQam();
        }

        /// <summary>
        /// Note: This does NOT include acquisition, tracking, AGC, or multipath deconvolution!
        /// 
        /// This shows how to...
        /// Create a QAM object with a specified constellation
        /// Optimize the constellation
        /// Visualize the constellation
        /// Convert data to states (upsample them to implement integration)
        /// Convert states to RF signal
        /// simulate transceiving by adding noise and convolution (multipath)
        /// (reception also reuires acquisition, tracking, deconvolution, and AGC, which must be done by app)
        /// Convert RF signal to states (downsample them to unimplement integration)
        /// Convert states to data
        /// Show BER stats
        /// 
        /// Integration is not implemented in the QAM class because it's easy to implement at the
        /// app level simply by upsampling the states before conversion to signal, and downsampling
        /// the states after conversion from signal.
        /// 
        /// </summary>
        /// <param name="noise"></param>
        /// <param name="numWordBits"></param>
        /// <param name="constellationTriangles"></param>
        public static void TestQam(
            double noise = 0.2,
            int numWordBits = 4,
            bool constellationTriangles = true)
        {
            var maxWordVal = (int)(Math.Pow(2, numWordBits) - .5);

            // create the QAM modulator/demodulator
            var qam = new QAM(QAM.QAM16_rectangular);

            // create a custom constellation
            qam.CreateConstellation(maxWordVal + 1, constellationTriangles, .5, .5);
            qam.Visualize();

            // optimize it
            Console.WriteLine("Optimizing symbol layout...");
            qam.OptimizeSymbolLayout(5);
            Console.WriteLine("Optimization complete");
            qam.Visualize();

            var totalBits = 0;
            var totalBitErrors = 0;

            while (true)
            {
                // Create a random packet
                var data = new List<int>();
                for (int n = 0; n < 8; n++)
                {
                    data.Add(Rnd.Next(0, maxWordVal + 1));
                    totalBits += numWordBits;
                }

                // convert data to QAM states
                var states = qam.SymbolsToStates(data);

                // convert states to an RF signal
                var signal = qam.StatesToSignal(states);

                // at this point we could upsample the signal to implement integration

                // simulate the radio transceiver (i.e. damage the signal)
                for (int n = 0; n < signal.Count; n++)
                {
                    signal[n] += NextGaussian() * noise;
                    // (Note that practical damage also includes some multipath, so we could convolute, here)
                }

                // note: you would need to acquire, track, and AGC the signal at this point

                // at this point we could downsample the signal to remove the integration

                // convert signal to states
                var rcvdStates = qam.SignalToStates(signal);

                // convert states to symbols
                var (rcvdData, avgDeviationSqr, maxDeviationSqr) = qam.StatesToSymbols(rcvdStates);

                // for every rcvd symbol
                for (int n = 0; n < rcvdData.Count; n++)
                {
                    if (data[n] != rcvdData[n])
                    {
                        // calc how many bits were in error
                        for (int b = 0; b < numWordBits; b++)
                        {
                            var bitMap = 1 << b;
                            if ((data[n] & bitMap) != (rcvdData[n] & bitMap))
                            {
                                totalBitErrors++;
                            }
                        }
                    }
                }
                Console.Write($"noise:{noise}  BER = {totalBitErrors / (double)totalBits}                               \r");
            }
        }

        public static double NextGaussian()
        {
            var u1 = Rnd.NextDouble();
            var u2 = Rnd.NextDouble();

            var rand_std_normal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                                Math.Sin(2.0 * Math.PI * u2);

            return rand_std_normal;
        }
    }
}
