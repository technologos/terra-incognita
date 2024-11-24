import carto

import logging
logger = logging.getLogger(__name__)

logging.basicConfig(filename='myapp.log', 
                    level=logging.INFO, 
                    format='%(asctime)s %(message)s')

random_seed = 0
pixel_width=600
pixel_height=600
pixels_per_tile=5

for random_seed in range(5):
    for grow_method in ['random', 'squared_inv_distances', 'squared_inv_centroid_distances']:
        logger.info(f'Started with random seed {random_seed}')
        test_map = carto.Map(pixel_width=pixel_width, pixel_height=pixel_height, pixels_per_tile=pixels_per_tile)
        test_map.setup(random_seed=random_seed)
        logger.info('Map setup complete')
        test_map.simulate_tectonics(n_plates=10, grow_method=grow_method)
        logger.info('Tectonics simulated')
        test_map.save(filename=f'test_map_{pixel_width}x{pixel_height}x{pixels_per_tile}_{random_seed}_{grow_method}.png', color_by='plate', show_starting_tile=True)
        logger.info(f'Finished with random seed {random_seed}')
