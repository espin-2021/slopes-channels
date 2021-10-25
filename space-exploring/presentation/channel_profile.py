# Instantiate main channel profiler
prf_main = ChannelProfiler(mg,
                           main_channel_only=True,
                           minimum_channel_threshold=3e4)# Instantiate multiple channel profiler
prf_multi = ChannelProfiler(mg,
                            main_channel_only=False,
                            minimum_channel_threshold=3e4)# Run channel profilers
prf_main.run_one_step()
prf_multi.run_one_step()# Plot channel profiles
prf_main.plot_profiles(ylabel="Elevation")
plt.show()prf_multi.plot_profiles(ylabel="Elevation")
plt.show()prf_main.plot_profiles_in_map_view()
plt.show()prf_multi.plot_profiles_in_map_view()
plt.show()
