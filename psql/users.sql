-- Query 1
SELECT * FROM nsrr_cloud_view

-- Query 2
SELECT user_id,user_id || '-' || authentication_token AS auth_token,dataset_slug 
FROM nsrr_cloud_view
ORDER BY user_id ASC;
