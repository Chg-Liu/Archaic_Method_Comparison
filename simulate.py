#!/usr/bin/env python3
import msprime
import tspop
import demesdraw
import numpy as np
import sys
import argparse
import pandas as pd


def initialize_demographic_model(N_ancestral=10000, N_homo=10000, N_chim=1000, N_AMH=10000, N_archaic=2800, N_AFR=14000, N_EUR_EAS=10000,
                                 T_Homo_Chim=230000,T_AMH_Archaic=24000, T_AFR_EUR_EAS=2000,
                                 T_introgression=1800, introgression_rate=0.05, introgression_duration=1,
                                 N_Bottleneck_EUR_EAS=None, T_Bottleneck_EUR_EAS=None,
                                 Bottleneck_EUR_EAS_Duration=None, T_EUR_EAS_growth=None, EUR_EAS_growth_rate=None):
    """
    Initialize a demographic model with the following parameters:
    :param N_ancestral: int, ancestral population size
    :param N_homo: int, homo population size
    :param N_chim: int, chimpanzee  population size
    :param N_archaic: int, archaic population size
    :param N_AMH: int, population size of ancestral modern human population
    :param N_AFR: int, African population size
    :param N_EUR_EAS: int, Eurasian population size
    :param T_Homo_Chim: int, time of split between homo and chimpanzee (number of generations ago)
    :param T_AMH_Archaic: int, time of split between archaic and modern humans (number of generations ago)
    :param T_AFR_EUR_EAS: int, time of out-of-Africa migration (number of generations ago)
    :param T_introgression: int, time of introgression (number of generations ago)
    :param introgression_rate: float, rate of introgression from Archaic in EUR_EAS
    :param introgression_duration: int, number of generations introgression lasts
    :param N_Bottleneck_EUR_EAS: int, Eurasian population size during bottleneck
    :param T_Bottleneck_EUR_EAS: int, time of the beginning of Eurasian bottleneck (number of generations ago)
    :param Bottleneck_EUR_EAS_Duration: int, duration of Eurasian bottleneck (number of generations)
    :param T_EUR_EAS_growth: int, time of Eurasian growth onset (number of generations ago)
    :param EUR_EAS_growth_rate: float, Eurasian growth rate
    :return: msprime.Demography, demographic model
    """
    demography = msprime.Demography()
    demography.add_population(name="Ancestral", initial_size=N_ancestral)
    demography.add_population(name="Homo", initial_size=N_homo)
    demography.add_population(name="Chim", initial_size=N_chim)
    demography.add_population(name="Archaic", initial_size=N_archaic)
    demography.add_population(name="AMH", initial_size=N_AMH)
    demography.add_population(name="AFR", initial_size=N_AFR)
    demography.add_population(name="EUR_EAS", initial_size=N_EUR_EAS)
    # start of introgression backward in time (i.e., end forward in time)
    demography.add_migration_rate_change(time=T_introgression - introgression_duration, source='EUR_EAS', dest="Archaic",
                                         rate=introgression_rate)
    # end of introgression backward in time (i.e., start forward in time)
    demography.add_migration_rate_change(time=T_introgression,
                                         source='EUR_EAS', dest="Archaic", rate=0)

    # TODO confirm that this is sufficient to reconstruct ground truth set of introgressed segments with tspop
    #  https://tspop.readthedocs.io/en/latest/simulationsetup.html
    demography.add_census(time=T_introgression + 0.01)
    # OOA
    demography.add_population_split(time=T_AFR_EUR_EAS, derived=["EUR_EAS", "AFR"], ancestral="AMH")
    # split of archaic and modern humans
    demography.add_population_split(time=T_AMH_Archaic, derived=["AMH", "Archaic"], ancestral="Homo")
    # split of archaic and modern humans
    demography.add_population_split(time=T_Homo_Chim, derived=["Homo", "Chim"], ancestral="Ancestral")

    # model bottleneck in EUR_EAS if specified
    if N_Bottleneck_EUR_EAS is not None:
        assert T_Bottleneck_EUR_EAS is not None and Bottleneck_EUR_EAS_Duration is not None, \
            "If N_bottleneck_eur_eas is set, T_bottleneck_eur_eas and bottleneck_eur_eas_duration must also be set."
        # start of bottleneck backward in time (i.e., end forward in time)
        demography.add_population_parameters_change(time=T_Bottleneck_EUR_EAS - Bottleneck_EUR_EAS_Duration,
                                                    population="EUR_EAS", initial_size=N_Bottleneck_EUR_EAS)
        # end of bottleneck backward in time (i.e., start forward in time)
        demography.add_population_parameters_change(time=T_Bottleneck_EUR_EAS, population="EUR_EAS",
                                                    initial_size=N_EUR_EAS)

    # model recent population growth in EUR_EAS if specified
    if T_EUR_EAS_growth is not None:
        assert EUR_EAS_growth_rate is not None, "If T_EUR_EAS_growth is set, EUR_EAS_growth_rate must also be set."
        # the population growths from N_EUR_EAS to it's
        # final size (N_EUR_EAS * np.exp(EUR_EAS_growth_rate * T_EUR_EAS_growth)) from T_EUR_EAS to 0 forward in time
        demography.add_population_parameters_change(time=0,
                                                    initial_size=N_EUR_EAS * np.exp(EUR_EAS_growth_rate *
                                                                                    T_EUR_EAS_growth),
                                                    growth_rate=EUR_EAS_growth_rate, population="EUR_EAS")
        # end of population growth backward in time (i.e., start forward in time)
        demography.add_population_parameters_change(time=T_EUR_EAS_growth, initial_size=N_EUR_EAS,
                                                    growth_rate=0, population="EUR_EAS")

    # make sure the events are in order
    demography.sort_events()
    return demography


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--N_ancestral", type=int, default=10000, help="Ancestral population size")
    parser.add_argument("--N_homo", type=int, default=10000, help="Homo population size")
    parser.add_argument("--N_chim", type=int, default=1000, help="Chimpanzee population size")
    parser.add_argument("--N_AMH", type=int, default=10000,
                        help="Population size of ancestral modern human population")
    parser.add_argument("--N_archaic", type=int, default=2800, help="Archaic population size")
    parser.add_argument("--N_AFR", type=int, default=14000, help="African population size")
    parser.add_argument("--N_EUR_EAS", type=int, default=10000, help="Eurasian population size")
    parser.add_argument("--T_Homo_Chim", type=int, default=230000,
                        help="Time of split between archaic and modern humans (number of generations ago)")
    parser.add_argument("--T_AMH_Archaic", type=int, default=24000,
                        help="Time of split between archaic and modern humans (number of generations ago)")
    parser.add_argument("--T_AFR_EUR_EAS", type=int, default=2000,
                        help="Time of out-of-Africa migration (number of generations ago)")
    parser.add_argument("--T_introgression", type=int, default=1800,
                        help="Time of introgression (number of generations ago)")
    parser.add_argument("--introgression_rate", type=float, default=0.05,
                        help="Rate of introgression from Archaic in EUR_EAS")
    parser.add_argument("--introgression_duration", type=int, default=1,
                        help="Number of generations introgression lasts")
    parser.add_argument("--N_Bottleneck_EUR_EAS", type=int, default=None,
                        help="Eurasian population size during bottleneck")
    parser.add_argument("--T_Bottleneck_EUR_EAS", type=int, default=None,
                        help="Time of the beginning of Eurasian bottleneck (number of generations ago)")
    parser.add_argument("--Bottleneck_EUR_EAS_Duration", type=int, default=None,
                        help="Duration of Eurasian bottleneck (number of generations)")
    parser.add_argument("--T_EUR_EAS_growth", type=int, default=None,
                        help="Time of Eurasian growth onset (number of generations ago)")
    parser.add_argument("--EUR_EAS_growth_rate", type=float, default=None,
                        help="Eurasian growth rate")
    parser.add_argument("--num_samples_afr", type=int, default=100,
                        help="Number of African samples to simulate")
    parser.add_argument("--num_samples_eur_eas", type=int, default=100,
                        help="Number of Eurasian samples to simulate")
    parser.add_argument("--num_samples_archaic", type=int, default=1,
                        help="Number of archaic samples to simulate")
    parser.add_argument("--num_samples_chim", type=int, default=1,
                        help="Number of chimpanzee samples to simulate")
    parser.add_argument("--archaic_sampling_time", type=int, default=3000,
                        help="Time of archaic sampling (number of generations ago)")
    parser.add_argument("--mutation_rate", type=float, default=1.2e-8,
                        help='Mutation rate per bp per generation')
    parser.add_argument("--recombination_rate", type=float, default=1e-8,
                        help='Recombination rate per bp per generation')
    parser.add_argument("--sequence_length", type=int, default=1e6,
                        help='Length of the simulated sequence')
    parser.add_argument("--output_trees", type=str, default="output.trees",
                        help='Output file for simulated treesequence')
    parser.add_argument("--output_vcf", type=str, default="output.vcf",
                        help='Output file for simulated VCF')
    parser.add_argument("--output_ancestry", type=str, default="output.ancestry_table",
                        help='Output file for population ancestry')
    parser.add_argument("--output_demes_figure", type=str, default=None,
                        help='If set plot demographic model using demesdraw and save to specified file.')
    parser.add_argument("--seed", type=int, default=42, help="Random seed for simulation")

    args = parser.parse_args()
    # initialize demographic model
    demography = initialize_demographic_model(args.N_ancestral, args.N_homo, args.N_chim, args.N_AMH, args.N_archaic, args.N_AFR, args.N_EUR_EAS,
                                              args.T_Homo_Chim, args.T_AMH_Archaic, args.T_AFR_EUR_EAS,args.T_introgression,
                                              args.introgression_rate, args.introgression_duration,
                                              args.N_Bottleneck_EUR_EAS, args.T_Bottleneck_EUR_EAS,
                                              args.Bottleneck_EUR_EAS_Duration, args.T_EUR_EAS_growth,
                                              args.EUR_EAS_growth_rate)

    if args.output_demes_figure is not None:
        graph = demography.to_demes()
        ax = demesdraw.tubes(graph)
        ax.figure.savefig(args.output_demes_figure)

    # define sample sizes and times
    sample_set = [msprime.SampleSet(args.num_samples_afr, population="AFR", time=0, ploidy=2),
                  msprime.SampleSet(args.num_samples_eur_eas, population="EUR_EAS", time=0, ploidy=2),
                  msprime.SampleSet(args.num_samples_archaic, population="Archaic", time=args.archaic_sampling_time,ploidy=2),
                  msprime.SampleSet(args.num_samples_chim, population="Chim", time=0, ploidy=2)]
    # do simulations
    # TODO allow setting a proper genetic map
    ts = msprime.sim_ancestry(samples=sample_set, demography=demography, random_seed=args.seed,
                              sequence_length=args.sequence_length, recombination_rate=args.recombination_rate)
    pa = tspop.get_pop_ancestry(ts, census_time=args.T_introgression+0.01)
    print(pa.ancestry_table)
    pa.ancestry_table.to_csv(args.output_ancestry,sep='\t',index=False)

    # add mutations
    mts = msprime.sim_mutations(ts, rate=args.mutation_rate, random_seed=args.seed)
    # save tree sequence
    mts.dump(args.output_trees)
    # write vcf
    with open(args.output_vcf, "w") as vcf_file:
        mts.write_vcf(vcf_file)
    vcf_file.close()


if __name__ == '__main__':
    main(sys.argv[1:])
