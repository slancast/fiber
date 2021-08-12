module load google-cloud-sdk
#!/usr/bin/env bash
for i in $(seq 103 154)
do
bash_script="humann2_gcloud_copy_$i"
chmod +x $bash_script
instance_name="instance-$i"
echo $instance_name
gcloud compute disks create $instance_name --source-snapshot humann2 --zone=us-central1-a
echo "disk created"
gcloud compute instances create $instance_name --custom-cpu=6 --custom-memory=39 --zone=us-central1-a --disk name=$instance_name,boot=yes 
echo "instance created"
gcloud compute scp $bash_script slancast@"$instance_name:~/" --zone=us-central1-a
bash_script2="\"./"$bash_script"\""
gcloud compute ssh slancast@$instance_name --zone=us-central1-a --command="screen -d -m bash -c $bash_script2" 
done

