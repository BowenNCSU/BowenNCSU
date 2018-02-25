import cProfile
import pstats

# Run profiler and save to twitter.stats
cProfile.run('import twitter', 'twitter.stats')
stats = pstats.Stats('twitter.stats')
# Clean up filenames for the report
stats.strip_dirs()
# Sort the statistics by the cumulative time spent
stats.sort_stats('cumulative')
stats.print_stats()
